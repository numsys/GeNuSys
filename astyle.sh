#!/bin/sh
TMP_PWD=$pwd
ASTYLE_OPTS="--style=allman --indent=spaces=4 --indent-classes --indent-switches --indent-cases --indent-namespaces --indent-labels --indent-preprocessor --indent-col1-comments --min-conditional-indent=2 --max-instatement-indent=80 --pad-oper --pad-header --unpad-paren --recursive --lineend=linux --unpad-paren --break-closing-brackets --convert-tabs --align-pointer=type --suffix=none --align-reference=type --indent-preproc-define --close-templates"
cd $1
astyle $ASTYLE_OPTS "*.c" $2 $3 $4 $5
astyle $ASTYLE_OPTS "*.hpp" $2 $3 $4 $5
astyle $ASTYLE_OPTS "*.cpp" $2 $3 $4 $5
astyle $ASTYLE_OPTS "*.h" $2 $3 $4 $5
astyle $ASTYLE_OPTS "*.ih" $2 $3 $4 $5
cd $TMP_PWD
