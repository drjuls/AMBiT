#! /bin/sh

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

echo "// Auto-generated source/build information."
echo "#define GIT_SOURCE_DESCRIPTION \"`echo_git_info`\""
echo "#define COMPILE_DATE \"`date`\""