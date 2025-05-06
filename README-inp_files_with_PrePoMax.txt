Seems that only C3D4 (tet elements) mesh works...

Bei inp-file parsing:
im inp-file müssen ein paar Dinge manuell geändert werden:
*Node -> *Node, NSET=Nall

Loads:
Leere Lasten löschen
*Cload
Internal_Selection-1_Concentrated_Force-1, 2, -50    ->    1, 2, -50 (CLOAD-Nodes müssen explizit eingetragen werden)