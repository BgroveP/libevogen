# Publish the staged artifact (artifact_staging/, produced by
# tools/rebuild_libs.sh and tools/build_wrapper.jl) as an in-repo
# artifact:
#
#   1. copy the staged content to artifacts/libevogen/,
#   2. git add it and read the staged git tree hash,
#   3. bind that hash in this package's Artifacts.toml,
#   4. install the same content into the local Julia depot.
#
# Usage:
#   julia tools/install_artifact.jl
#
# Afterwards commit artifacts/libevogen/ and Artifacts.toml. Pkg resolves
# an artifact whose git-tree-sha1 matches a tree inside the package's
# own repository directly from that repository, so no download URL is
# needed; this works for private and public repositories alike.
using Pkg.Artifacts: bind_artifact!, artifact_path
import Pkg.Artifacts
import Base.BinaryPlatforms: Platform

const PKG_ROOT = dirname(@__DIR__)
const STAGING = joinpath(PKG_ROOT, "artifact_staging")
const IN_REPO_ARTIFACT = joinpath(PKG_ROOT, "artifacts", "libevogen")
const ARTIFACTS_TOML = joinpath(PKG_ROOT, "Artifacts.toml")

isdir(joinpath(STAGING, "lib")) || error(
    "no staged artifact found; run tools/rebuild_libs.sh first")
isfile(joinpath(STAGING, "lib", "libevogen_wrapper.so")) || error(
    "the staged artifact has no wrapper; run tools/build_wrapper.jl first")

@info "Publishing in-repo artifact" dst = IN_REPO_ARTIFACT
rm(IN_REPO_ARTIFACT; recursive = true, force = true)
mkpath(dirname(IN_REPO_ARTIFACT))
cp(STAGING, IN_REPO_ARTIFACT)

# Stage the content in git and take the tree hash of artifacts/libevogen
# from the index (git write-tree --prefix prints the subtree hash of a
# staged directory). This is the hash that Pkg must find in the
# repository, so the files must be committed exactly as added (beware of
# autocrlf rewriting text files - the headers here carry their original
# CRLF endings and .gitattributes is not used).
run(`git add artifacts/libevogen`)
hash = Base.SHA1(strip(read(`git write-tree --prefix=artifacts/libevogen/`, String)))
@info "Staged tree hash" hash = string(hash)

# Bind it in Artifacts.toml for linux x86_64 (glibc).
platform = Platform("x86_64", "linux"; libc = "glibc")
bind_artifact!(ARTIFACTS_TOML, "libevogen", hash;
               platform = platform, lazy = false, force = true)
@info "Artifact bound in Artifacts.toml" artifacts_toml = ARTIFACTS_TOML

# Install the same content into the local depot so the package can be
# used (and tested) before the next release.
dest = artifact_path(hash)
rm(dest; recursive = true, force = true)
mkpath(dirname(dest))
cp(IN_REPO_ARTIFACT, dest)
@info "Artifact installed" location = dest

println("""

Commit the artifact and the manifest:

    git add Artifacts.toml
    git commit -m "Rebuild libevogen artifact ($(string(hash)[1:8]))"

Consumers get the new artifact on their next package update.
""")
