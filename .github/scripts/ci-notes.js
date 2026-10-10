#!/usr/bin/env node
/**
 * 生成「上一个 Release → 本次构建」之间的提交清单，供 gh release create --notes
 * 前置到 GitHub 自动生成的发布说明之前。
 *
 * 为什么需要它：GitHub 自动生成的发布说明只统计「已合并的 Pull Request」，
 * 直接推送到分支的提交不会出现在里面。本脚本补上这段区间内的全部提交。
 *
 * 用法: node .github/scripts/ci-notes.js <本次将要创建的 tag>
 * 依赖: 环境变量 GITHUB_REPOSITORY、GITHUB_SHA；GITHUB_TOKEN 可选（用于提高配额）。
 * 输出: Markdown 片段写到 stdout；找不到基线或发生异常时输出为空，不阻断发布。
 *
 * 注意: 本次要创建的 tag 此时还不存在，所以比对区间的终点用的是 GITHUB_SHA
 *       所指向的提交，而不是 tag 名。
 */

const API = 'https://api.github.com';
const repo = process.env.GITHUB_REPOSITORY || '';
const token = process.env.GITHUB_TOKEN || '';
const head = (process.env.GITHUB_SHA || '').trim();
const currentTag = (process.argv[2] || '').trim();

async function api(path) {
  const headers = {
    Accept: 'application/vnd.github+json',
    'X-GitHub-Api-Version': '2022-11-28',
    'User-Agent': 'molecular-viewer-ci',
  };
  if (token) {
    headers.Authorization = `Bearer ${token}`;
  }
  const res = await fetch(`${API}${path}`, { headers });
  if (!res.ok) {
    throw new Error(`GET ${path} -> ${res.status} ${res.statusText}`);
  }
  return res.json();
}

/**
 * 把 compare 接口返回的提交转成 Markdown 列表。
 * 跳过 merge commit —— 它们的内容已由各自被合并的提交体现。
 */
function formatNotes(baseTag, commits) {
  const lines = commits
    .filter((c) => (c.parents || []).length < 2)
    .map((c) => {
      const sha = c.sha.slice(0, 7);
      const subject = String(c.commit.message).split('\n')[0].trim();
      const who =
        c.author && c.author.login
          ? `@${c.author.login}`
          : String(c.commit.author.name || '').trim();
      return `* \`${sha}\` ${subject} — ${who}`;
    });

  if (lines.length === 0) {
    return '';
  }
  return [`## 提交记录（\`${baseTag}\` 之后）`, '', ...lines].join('\n');
}

async function main() {
  if (!repo || !head) {
    return '';
  }

  const releases = await api(`/repos/${repo}/releases?per_page=30`);
  // 最近一个既有 Release 作为基线；排除本次要创建的那个 tag（重跑场景）
  const baseTag = releases.map((r) => r.tag_name).find((t) => t && t !== currentTag);
  if (!baseTag) {
    return '';
  }

  const compare = await api(`/repos/${repo}/compare/${baseTag}...${head}`);
  return formatNotes(baseTag, compare.commits || []);
}

main()
  .then((notes) => {
    if (notes) {
      process.stdout.write(`${notes}\n`);
    }
  })
  .catch((err) => {
    // 提交清单是锦上添花，任何失败都不应该挡住发布
    console.error(`ci-notes: 跳过提交清单（${err.message}）`);
  });
