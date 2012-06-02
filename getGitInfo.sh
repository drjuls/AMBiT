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

echo "// Auto-generated source/build information."
echo "#define GIT_LAST_TAG \"`echo_git_tag`\""
echo "#define GIT_SOURCE_DESCRIPTION \"`echo_git_info`\""
echo "#define COMPILE_DATE \"`date`\""
