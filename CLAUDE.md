# umi-tools-rs

## Toolchain

All tools are managed by [mise](https://mise.jdx.dev/) — see `mise.toml` for versions.

**Important:** Claude Code must be launched from a shell where mise is activated (`mise activate` in shell profile). The Bash tool inherits the parent shell's environment snapshot. If `cargo` / `uv` are not found, restart Claude Code from a properly activated shell rather than using `mise exec` workarounds.

## Common commands

```sh
mise run test     # cargo test --workspace
mise run compat   # umi-tools compatibility tests (pytest)
mise run lint     # cargo clippy (strict)
mise run fmt      # cargo fmt --check
mise run ci       # all of the above
```
