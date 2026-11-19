#!/usr/bin/env bash
# Replace the 0.0.0-dev placeholder with the version of the pushed tag.
set -euo pipefail

VERSION="${GITHUB_REF_NAME#v}"
sed -i "s/^version = \"0.0.0-dev\"/version = \"$VERSION\"/" Cargo.toml crates/umi-core/Cargo.toml
sed -i "s|umi-core = { path = \"crates/umi-core\", version = \"0.0.0-dev\" }|umi-core = { path = \"crates/umi-core\", version = \"$VERSION\" }|" Cargo.toml
for manifest in Cargo.toml crates/umi-core/Cargo.toml; do
    grep -q "^version = \"$VERSION\"" "$manifest"
done
grep -q "version = \"$VERSION\" }" Cargo.toml
