#!/bin/sh

set -e

if ! type git > /dev/null 2>&1; then
    echo "ERROR: git not found, can't continue" >&2
    exit 1
fi

echo "Checking for tabs" >&2
! git --no-pager grep -InP --heading "\t" -- . ':!third_party/**/*' || ret=1

echo "Checking for carriage returns" >&2
! git --no-pager grep -InP --heading "\r" -- . ':!third_party/**/*' || ret=1

echo "Checking for trailing spaces" >&2
! git --no-pager grep -InP --heading " $" -- . ':!third_party/**/*' ':!*.patch' || ret=1

# Test only "new" commits, that is, commits that are not upstream on
# the default branch.
if git fetch -q https://gitlab.com/AOMediaCodec/SVT-AV1.git HEAD; then
    FETCH_HEAD=FETCH_HEAD
else
    # in case the fetch failed, maybe internet issues, try to resolve a local default branch's ref, checked-out or not
    FETCH_HEAD=$(git rev-parse refs/remotes/origin/HEAD)
fi

# default to master if we have no origin remote
: "${FETCH_HEAD:=master}"

if ! git merge-tree "$(git merge-base HEAD "$FETCH_HEAD")" HEAD "$FETCH_HEAD"; then
    echo "ERROR: failed to simulate a merge, check to see if a merge is possible" >&2
fi

if git diff --exit-code --diff-filter=d --name-only "^$FETCH_HEAD" > /dev/null 2>&1; then
    echo "No differences to upstream's default, skipping further tests"
    exit 0
fi

while read -r file; do
    if ! test -f "$file"; then
        printf "Ignoring file not found: '%s'\n" "$file"
        continue
    fi
    if test -n "$(tail -c1 "$file")"; then
        printf "No newline at end of %s\n" "$file"
        ret=1
    fi
done << EOF
$(
    git diff --name-only --diff-filter=d "$FETCH_HEAD" -- . \
        ':!third_party' ':!test/e2e_test/test_vector_list.txt' \
        ':!test/vectors/smoking_test.cfg' \
        ':!test/vectors/video_src.cfg' \
        ':!*.png'
)
EOF

while read -r i; do
    printf "Checking commit message of %s\n" "$i" >&2
    msg=$(git log --format=%B -n 1 "$i")
    if test -n "$(printf '%s' "$msg" | sed -n 2p)"; then
        printf "Malformed commit message in %s, second line must be empty\n" "$i"
        ret=1
    fi
    if printf '%s' "$msg" | head -1 | grep -q '\.$'; then
        printf "Malformed commit message in %s, trailing period in subject line\n" "$i"
        ret=1
    fi
done << EOF
$(git rev-list HEAD "^$FETCH_HEAD")
EOF
exit ${ret:-0}
