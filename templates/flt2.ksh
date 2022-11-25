#!/bin/ksh
#-*-ksh-*-
editfst -s $1 -d $2 -l << EOF
exclure(-1,['TM', 'UW', 'VW', 'UUW','VVW'], -1, -1, 26314400, -1, -1) 
desire(-1,['^^'])
desire(-1,['>>'])
desire(-1,['^>'])
desire(-1,['TM','UW', 'VW', 'UUW','VVW'],-1,-1,-1,-1,-1)
EOF
