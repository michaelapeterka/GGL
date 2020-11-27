rule[
 ruleID "r3: imine-enamine tautomerization"
 left[
  node[id 0 label "H"]
  node[id 4 label "H+"]
  edge[source 0 target 1 label "-"]
  edge[source 1 target 2 label "-"]
  edge[source 2 target 3 label "="]
 ]
 context[
  node[id 1 label "C"]
  node[id 2 label "C"]
  node[id 3 label "N"]
  node[id 5 label "*"]
  edge[source 1 target 5 label "-"]
 ]
 right[
  node[id 0 label "H+"]
  node[id 4 label "H"]
  edge[source 1 target 2 label "="]
  edge[source 2 target 3 label "-"]
  edge[source 3 target 4 label "-"]
 ]
] 
