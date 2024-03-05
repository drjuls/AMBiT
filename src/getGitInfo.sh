#! /bin/bash

function echo_git_info {
    if [[ $(which git 2> /dev/null) ]]
    then
        local STATUS
        STATUS=$(git status 2>/dev/null)
        if [[ -z $STATUS ]]
        then
            return
        fi
        echo "`git symbolic-ref HEAD 2> /dev/null | cut -b 12-` `git log --pretty=format:\"%h\" -1`"
        return
    fi
}

function echo_git_tag {
    if [[ $(which git 2> /dev/null) ]]
    then
        local STATUS
        STATUS=$(git status 2>/dev/null)
        if [[ -z $STATUS ]]
        then
            return
        fi
        echo "`git describe --tags --abbrev=0 HEAD 2> /dev/null`"
        return
    fi
}

function echo_git_diff {
    if [[ $(which git 2> /dev/null) ]]
    then
        local STATUS
        STATUS=$(git status 2>/dev/null)
        if [[ -z $STATUS ]]
        then
            return
        fi
#
# the following sed commands mean, in order:
# 1. Replace '\' with '\\'
# 2. Add '\n \' to end of every file
# 3. Remove the final '\' from the last line
# 4. Replace " with \"
#
        echo "`git diff -w --abbrev=0 HEAD 2> /dev/null | sed -e 's/\\\\/\\\\\\\\/g' -e 's/$/\\\\n\\\\/' -e '$s/\\\\$//' -e 's/\"/\\\\"/g' `"
        return
    fi
}

echo "// Auto-generated source/build information."
echo "#define GIT_LAST_TAG \"`echo_git_tag`\""
echo "#define GIT_SOURCE_DESCRIPTION \"`echo_git_info`\""
echo "#define COMPILE_DATE \"`date`\""
echo "#define GIT_DIFF_OUTPUT \"`echo_git_diff`\""
