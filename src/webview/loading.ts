// Loading helpers shared by the provider and the tool paths: text decoding,
// chunked reads and the determinate progress notification.
import * as vscode from 'vscode';
import * as fs from 'fs';

/**
 * Decode a text file the way the writers on other machines wrote it.
 *
 * The readers here assume UTF-8, but a Windows job editor (and some
 * crystallographic tools) save `.gjf`, `.inp` and `.data` as UTF-16LE with a
 * BOM. Decoding those bytes as UTF-8 turns every character into a replacement
 * or a NUL and no parser can read the file. The BOM decides the encoding;
 * everything else keeps the previous UTF-8 default (TextDecoder strips a UTF-8
 * BOM itself).
 */
export function decodeText(bytes: Uint8Array): string {
    if (bytes.length >= 2 && bytes[0] === 0xff && bytes[1] === 0xfe) {
        return new TextDecoder('utf-16le').decode(bytes);
    }
    return new TextDecoder().decode(bytes);
}

/** The progress bar a file load uses while it reads and parses. */
export type LoadProgress = vscode.Progress<{ message?: string; increment?: number }>;

/**
 * A file this large is slow enough to read and parse that its own progress
 * notification is worth showing. Smaller structures open instantly, and a
 * notification that flashes for every one of them is worse than none.
 */
export const LOAD_PROGRESS_MIN_BYTES = 4 * 1024 * 1024;
/**
 * A log is watched at a much smaller size: its single pass builds a frame per
 * optimization step, and a 1 MB optimisation log is already a visible wait.
 * The threshold only exists to keep a 5 KB output from flashing a notification.
 */
export const LOAD_PROGRESS_MIN_LOG_BYTES = 256 * 1024;
/** Bytes per read, so a cancel is felt and the percentage can move. */
export const LOAD_PROGRESS_READ_CHUNK = 1 << 20;
/** Atoms a page payload needs before building it is slow enough to report. */
export const VIEW_PROGRESS_MIN_ATOMS = 250000;
/** A reporter for the paths that chose not to open a notification. */
export const NO_PROGRESS: LoadProgress = { report: () => { /* not shown */ } };

export function formatBytes(bytes: number): string {
    if (bytes < 1024) return bytes + ' B';
    if (bytes < 1024 * 1024) return (bytes / 1024).toFixed(1) + ' KB';
    if (bytes < 1024 * 1024 * 1024) return (bytes / (1024 * 1024)).toFixed(1) + ' MB';
    return (bytes / (1024 * 1024 * 1024)).toFixed(2) + ' GB';
}

/** A file's size, or 0 when the file system cannot state one. */
export async function fileSize(uri: vscode.Uri): Promise<number> {
    try {
        return (await vscode.workspace.fs.stat(uri)).size || 0;
    } catch {
        return 0;
    }
}

/**
 * Whether opening this file is expected to take long enough to report.
 *
 * A log is checked at half the size because its single pass costs more per
 * byte. A VMD script is not decided here: the script itself is tiny, so its
 * own branch measures the structures it names before loading them.
 */
export function wantsLoadProgress(size: number, ext: string): boolean {
    if (ext === 'log' || ext === 'out') return size >= LOAD_PROGRESS_MIN_LOG_BYTES;
    return size >= LOAD_PROGRESS_MIN_BYTES;
}

/** One turn of the host message loop, so a progress message is painted first. */
export function yieldToHost(): Promise<void> {
    return new Promise<void>(resolve => setTimeout(resolve, 0));
}

/**
 * Map a 0..1 fraction onto the [from, to] slice of one progress bar.
 *
 * The workbench draws a determinate bar from integer percentages: a fractional
 * increment such as 0.3 is invisible, and a first report of 0 can leave the
 * notification indeterminate. Each step is therefore a whole percentage point,
 * and the first report of a slice carries one point of its own so the bar is
 * determinate from the start. A phase label with no point left to spend is
 * still sent as a message; the bar is determinate by then.
 */
export function progressRange(
    progress: LoadProgress,
    from: number,
    to: number
): (fraction: number, message?: string) => void {
    const start = Math.round(from);
    const end = Math.round(to);
    let last = start;
    let announced = false;
    let lastMessage: string | undefined;
    const emit = (increment: number, message?: string): void => {
        if (increment === 0) {
            // A label that has not changed is not news; repeating it once per
            // parser line would flood the notification for no visible gain.
            if (message === undefined || message === lastMessage) return;
            progress.report({ message });
        } else {
            progress.report(message === undefined ? { increment } : { increment, message });
        }
        if (message !== undefined) lastMessage = message;
    };
    return (fraction, message) => {
        const clamped = Math.max(0, Math.min(1, isFinite(fraction) ? fraction : 0));
        const target = Math.min(end, start + Math.floor((end - start) * clamped));
        if (!announced && message !== undefined) {
            announced = true;
            if (last < end) {
                last += 1;
                emit(1, message);
            } else {
                emit(0, message);
            }
            if (target <= last) return;
        }
        if (target <= last) {
            if (message !== undefined && target === last) emit(0, message);
            return;
        }
        emit(target - last, message);
        last = target;
    };
}

/**
 * Run a long load under a cancellable progress notification.
 *
 * The editor's own cancellation token and the notification's Cancel button are
 * folded into one flag: closing the tab and pressing Cancel both stop the load
 * at its next check. A synchronous parse cannot be interrupted line by line, so
 * the result is discarded at the next check rather than written into a document
 * nobody is waiting for.
 */
export async function withLoadProgress<T>(
    title: string,
    token: vscode.CancellationToken | undefined,
    work: (progress: LoadProgress, isCancelled: () => boolean) => Promise<T>
): Promise<T> {
    let cancelled = !!token && token.isCancellationRequested;
    return vscode.window.withProgress({
        location: vscode.ProgressLocation.Notification,
        title,
        cancellable: true
    }, async (progress, progressToken) => {
        const subscriptions: vscode.Disposable[] = [];
        if (token) {
            subscriptions.push(token.onCancellationRequested(() => { cancelled = true; }));
        }
        subscriptions.push(progressToken.onCancellationRequested(() => { cancelled = true; }));
        // A cancelled notification no longer accepts reports; the work may still
        // be finishing its current pass, and a late report must not throw.
        const guarded: LoadProgress = {
            report: value => {
                if (!cancelled) progress.report(value);
            }
        };
        try {
            const result = await work(guarded, () => cancelled);
            if (cancelled) throw new vscode.CancellationError();
            return result;
        } finally {
            for (const subscription of subscriptions) subscription.dispose();
        }
    });
}

/**
 * Run one heavy parse (or read) behind the same notification when the file is
 * large enough; a small file takes the plain path with no notification.
 */
export function withToolProgress<T>(
    title: string,
    heavy: boolean,
    work: (progress: LoadProgress, isCancelled: () => boolean) => Promise<T>
): Promise<T> {
    return heavy ? withLoadProgress(title, undefined, work) : work(NO_PROGRESS, () => false);
}

/**
 * Read a local file in chunks, reporting how much of it is in.
 *
 * `workspace.fs.readFile` returns one buffer with no interim feedback, so a
 * 500 MB log sits there for seconds with nothing to show. The bytes still end
 * up in one buffer — that is what the parsers need — but they are read in
 * chunks so the notification can name a percentage and a Cancel is felt
 * between them. Anything that is not a local file keeps the workspace API: a
 * virtual scheme has no path `fs` could open.
 */
export async function readFileWithProgress(
    uri: vscode.Uri,
    size: number,
    report: (fraction: number, message: string) => void,
    isCancelled: () => boolean
): Promise<Uint8Array> {
    report(0, size > 0 ? 'Reading ' + formatBytes(size) + ' · 0%' : 'Reading the file');
    if (uri.scheme !== 'file' || size <= 0) {
        const bytes = await vscode.workspace.fs.readFile(uri);
        if (isCancelled()) throw new vscode.CancellationError();
        report(1, 'Reading ' + formatBytes(bytes.length) + ' · 100%');
        return bytes;
    }
    const handle = await fs.promises.open(uri.fsPath, 'r');
    try {
        let capacity = Math.max(size, 1);
        let buffer = Buffer.allocUnsafe(capacity);
        let read = 0;
        let lastPercent = -1;
        for (;;) {
            if (isCancelled()) throw new vscode.CancellationError();
            if (read === capacity) {
                // The file grew after its stat: probe one byte before paying for
                // a second buffer, so an exactly-sized file allocates only once.
                const probe = Buffer.allocUnsafe(1);
                const { bytesRead } = await handle.read(probe, 0, 1, read);
                if (bytesRead <= 0) break;
                const grown = Buffer.allocUnsafe(Math.max(capacity * 2, read + LOAD_PROGRESS_READ_CHUNK));
                buffer.copy(grown, 0, 0, read);
                grown[read] = probe[0];
                read += 1;
                buffer = grown;
                capacity = grown.length;
            }
            const { bytesRead } = await handle.read(
                buffer, read, Math.min(capacity - read, LOAD_PROGRESS_READ_CHUNK), read
            );
            if (bytesRead <= 0) break;
            read += bytesRead;
            const percent = Math.min(99, Math.floor((read / size) * 100));
            // One message per whole percent, not per chunk: a 1 GB file is 100
            // updates in total, and a 52 MB one still fills its bar step by step.
            if (percent !== lastPercent) {
                report(percent / 100, 'Reading ' + formatBytes(size) + ' · ' + percent + '%');
                lastPercent = percent;
            }
        }
        if (isCancelled()) throw new vscode.CancellationError();
        report(1, 'Reading ' + formatBytes(read) + ' · 100%');
        return read === buffer.length ? buffer : buffer.subarray(0, read);
    } finally {
        await handle.close();
    }
}
