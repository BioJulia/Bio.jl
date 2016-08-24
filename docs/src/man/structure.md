# Bio.Structure: Macromolecular Structures

```@meta
CurrentModule = Bio.Structure
DocTestSetup = quote
    using Bio.Structure
    path = Pkg.dir("Bio", "test", "BioFmtSpecimens")
    if !isdir(path)
        run(`git clone --depth 1 https://github.com/BioJulia/BioFmtSpecimens.git $(path)`)
    end
    filepath_1EN2 = Pkg.dir("Bio", "test", "BioFmtSpecimens", "PDB", "1EN2.pdb")
end
```

The `Bio.Structure` module provides functionality to read [Protein Data Bank](http://www.rcsb.org/pdb/home/home.do) (PDB) files and manipulate macromolecular structures. It is designed to be used as a platform on which others can build to create tools for macromolecular structure analysis.


## Parsing PDB files

To download a PDB file:

```julia
downloadpdb("1EN2")
```

To parse a PDB file into a Structure-Model-Chain-Residue-Atom framework:

```julia
julia> struc = read(filepath_1EN2, PDB)
Name                        -  1EN2.pdb
Number of models            -  1
Chain(s)                    -  A
Number of residues          -  85
Number of point mutations   -  5
Number of other molecules   -  5
Number of water molecules   -  76
Number of atoms             -  614
Number of hydrogens         -  0
Number of disordered atoms  -  27
```

The elements of `struc` can be accessed as follows:

| Command                     | Returns                                                                         | Return type      |
| :-------------------------- | :------------------------------------------------------------------------------ | :--------------- |
| `struc[1]`                  | Model 1                                                                         | `Model`          |
| `struc[1]['A']`             | Model 1, chain A                                                                | `Chain`          |
| `struc['A']`                | The lowest model (model 1), chain A                                             | `Chain`          |
| `struc['A']["50"]`          | Model 1, chain A, residue 50                                                    | `Residue`        |
| `struc['A'][50]`            | Shortcut to above if it is not a hetero residue and the insertion code is blank | `Residue`        |
| `struc['A']["H_90"]`        | Model 1, chain A, hetero residue 90                                             | `Residue`        |
| `struc['A'][50]["CA"]`      | Model 1, chain A, residue 50, atom name CA                                      | `Atom`           |
| `struc['A'][15]["CG"]['A']` | For disordered atoms, access a specific location                                | `DisorderedAtom` |
| `struc['A'][15]["CG"]`      | For disordered atoms, access the default location                               | `Atom`           |

Disordered atoms are stored in a `DisorderedAtom` container but calls fall back to the default atom, so disorder can be ignored if you are not interested in it.

Disordered residues (i.e. point mutations with different residue names) are stored in a `DisorderedResidue` container.

The idea is that disorder will only bother you if you want it to. See the [Biopython discussion](http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ#How_is_disorder_handled.3F) for more.

Properties can be retrieved as follows:

| Command                  | Returns                                          | Return type            |
| :----------------------- | :----------------------------------------------- | :--------------------- |
| `structurename(struc)`   | Name of a structure                              | `ASCIIString`          |
| `modelnumbers(struc)`    | Model numbers in a structure                     | `Array{Int64,1}`       |
| `chainids(struc)`        | Chain IDs in a structure                         | `Array{Char,1}`        |
| `modelnumber(model)`     | Number of a model                                | `Int64`                |
| `chainids(model)`        | Chain IDs in a model                             | `Array{Char,1}`        |
| `chainid(chain)`         | Chain ID of a chain                              | `Char`                 |
| `resids(chain)`          | Residue IDs in a chain                           | `Array{ASCIIString,1}` |
| `resname(res)`           | Name of a residue                                | `ASCIIString`          |
| `chainid(res)`           | Chain ID of a residue                            | `Char`                 |
| `resnumber(res)`         | Residue number of a residue                      | `Int64`                |
| `inscode(res)`           | Insertion code of a residue                      | `Char`                 |
| `ishetres(res)`          | `true` if the residue consists of hetero atoms   | `Bool`                 |
| `atomnames(res)`         | Atom names in a residue                          | `Array{ASCIIString,1}` |
| `isdisorderedres(res)`   | `true` if the residue has multiple residue names | `Bool`                 |
| `resid(res)`             | Residue ID of a residue                          | `ASCIIString`          |
| `resid(res; full=true)`  | Residue ID of a residue including chain          | `ASCIIString`          |
| `ishetatom(atom)`        | `true` if the atom is a hetero atom              | `Bool`                 |
| `serial(atom)`           | Serial number of an atom                         | `Int64`                |
| `atomname(atom)`         | Name of an atom                                  | `ASCIIString`          |
| `altlocid(atom)`         | Alternative location ID of an atom               | `Char`                 |
| `resname(atom)`          | Residue name of an atom                          | `ASCIIString`          |
| `chainid(atom)`          | Chain ID of an atom                              | `Char`                 |
| `resnumber(atom)`        | Residue number of an atom                        | `Int64`                |
| `inscode(atom)`          | Insertion code of an atom                        | `Char`                 |
| `x(atom)`                | x coordinate of an atom                          | `Float64`              |
| `y(atom)`                | y coordinate of an atom                          | `Float64`              |
| `z(atom)`                | z coordinate of an atom                          | `Float64`              |
| `coords(atom)`           | coordinates of an atom                           | `Array{Float64,1}`     |
| `occupancy(atom)`        | Occupancy of an atom (default is `1.0`)          | `Float64`              |
| `tempfac(atom)`          | Temperature factor of an atom (default is `0.0`) | `Float64`              |
| `element(atom)`          | Element of an atom (default is `""`)             | `ASCIIString`          |
| `charge(atom)`           | Charge of an atom (default is `""`)              | `ASCIIString`          |
| `isdisorderedatom(atom)` | `true` if the atom is disordered                 | `Bool`                 |
| `resid(atom)`            | Residue ID of an atom                            | `ASCIIString`          |
| `resid(atom; full=true)` | Residue ID of an atom including chain            | `ASCIIString`          |

The coordinates of an atom can be set using `x!`, `y!`, `z!` and `coords!`.


## Manipulating structures

Elements can be looped over to reveal the sub-elements in the correct order:

```julia
for model in struc
    for chain in model
        for res in chain
            for atom in res
                # Do something
            end
        end
    end
end
```

Models are ordered numerically; chains are ordered by character, except the empty chain is last; residues are ordered by residue number and insertion code with hetero residues after standard residues; atoms are ordered by atom serial.

`collect`, `collectresidues` and `collectatoms` can be used to get lists of sub-elements.

Selectors are functions passed as additional arguments to `collectresidues` and `collectatoms`. Only residues/atoms that return `true` when passed to the selector are retained.

| Command                                                 | Action                                                            | Return type                |
| :------------------------------------------------------ | :---------------------------------------------------------------- | :------------------------- |
| `collect(struc['A'][50])`                               | Collect the sub-elements of an element, e.g. atoms from a residue | `Array{AbstractAtom,1}`    |
| `collectresidues(struc)`                                | Collect the residues of an element                                | `Array{AbstractResidue,1}` |
| `collectatoms(struc) `                                  | Collect the atoms of an element                                   | `Array{AbstractAtom,1}`    |
| `collectatoms(struc, calphaselector)`                   | Collect the C-alpha atoms of an element                           | `Array{AbstractAtom,1}`    |
| `collectatoms(struc, calphaselector, disorderselector)` | Collect the disordered C-alpha atoms of an element                | `Array{AbstractAtom,1}`    |

It is easy to define your own selector. The below will collect all atoms with x coordinate less than 0:

```julia
xselector(atom::AbstractAtom) = x(atom) < 0
collectatoms(struc, xselector)
```

Alternatively, you can use an anonymous function:

```julia
collectatoms(struc, atom -> x(atom) < 0)
```

Selectors can also be passed to many other functions, including `read`:

```julia
julia> calphas = read(filepath_1EN2, PDB, calphaselector)
Name                        -  1EN2.pdb
Number of models            -  1
Chain(s)                    -  A
Number of residues          -  85
Number of point mutations   -  5
Number of other molecules   -  0
Number of water molecules   -  0
Number of atoms             -  85
Number of hydrogens         -  0
Number of disordered atoms  -  2
```

`countmodels`, `countchains`, `countresidues` and `countatoms` can be used to count elements. For example:

```@meta
DocTestSetup = quote
    using Bio.Structure
    filepath_1EN2 = Pkg.dir("Bio", "test", "BioFmtSpecimens", "PDB", "1EN2.pdb")
    struc = read(filepath_1EN2, PDB)
end
```

```julia
julia> countatoms(struc)
754

julia> countatoms(struc, calphaselector)
85

julia> countresidues(struc, stdresselector)
85
```

`organise`, `organisemodel` and `organisestructure` can be used to organise sub-elements into elements:

| Command                                  | Action                                    | Return type                |
| :--------------------------------------- | :---------------------------------------- | :------------------------- |
| `organise(collectatoms(struc))`          | Organise an atom list into a residue list | `Array{AbstractResidue,1}` |
| `organise(collectresidues(struc))`       | Organise a residue list into a chain list | `Array{Chain,1}`           |
| `organise(struc['A'])`                   | Organise chain(s) into a model            | `Model`                    |
| `organise(struc[1])`                     | Organise model(s) into a structure        | `ProteinStructure`         |
| `organisemodel(collectatoms(struc))`     | Organise elements into a model            | `Model`                    |
| `organisestructure(collectatoms(struc))` | Organise elements into a structure        | `ProteinStructure`         |

The sequence of a protein can be retrieved:

```julia
julia> AminoAcidSequence(struc, stdresselector)
85aa Amino Acid Sequence:
RCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENKCWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYRC
```


## Spatial calculations

Distances can be calculated. The minimum distance between residue 10 and 20 is:

```julia
julia> distance(struc['A'][10], struc['A'][20])
10.782158874733762
```

Bond angles and dihedral angles can be calculated:

```julia
julia> rad2deg(bondangle(struc['A'][50]["N"], struc['A'][50]["CA"], struc['A'][50]["C"]))
69.22234153916602

julia> rad2deg(dihedralangle(struc['A'][50]["N"], struc['A'][50]["CA"], struc['A'][50]["C"], struc['A'][51]["N"]))
-177.38288114072924

julia> rad2deg(psiangle(struc['A'][50], struc['A'][51]))
-177.38288114072924
```

Further spatial functions including `contactmap`, `ramachandranangles`, `rmsd` and `displacements` are described in the Examples section below.


## Writing PDB files

PDB format files can be written:

```julia
writepdb("1EN2_out.pdb", struc)
```

Any element type can be given as input. Selectors can also be given as additional arguments:

```julia
writepdb("1EN2_out.pdb", struc, backboneselector)
```


## Examples

A few examples of `Bio.Structure` usage are given below.

**A)** To plot the temperature factors of a protein, if you have Gadfly installed:

```julia
using Gadfly
calpha_atoms = collectatoms(struc, calphaselector)
res_numbers = map(resnumber, calpha_atoms)
temp_facs = map(tempfac, calpha_atoms)
plot(x=res_numbers,
    y=temp_facs,
    Guide.xlabel("Residue number"),
    Guide.ylabel("Temperature factor"),
    Geom.line)
```

**B)** To find all C-alpha atoms within 5 Angstroms of residue 38:

```julia
for atom in calpha_atoms
    if distance(struc['A'][38], atom) < 5.0 && resnumber(atom) != 38
        show(atom)
    end
end
```

**D)** To show the contact map of a structure, if you have Hinton installed:

```julia
using Hinton
cbeta_atoms = collectatoms(struc, cbetaselector)
contacts = contactmap(cbeta_atoms, 7.0)
println(hintontxt(contacts))
```

`cbetaselector` selects C-beta atoms, or C-alpha atoms for glycine residues. `contactmap` can also be given two structural elements as arguments, in which case a non-symmetrical 2D array is returned showing contacts between the elements.

**E)** To show the Ramachandran phi/psi angle plot of a structure, if you have Gadfly installed:

```julia
using Gadfly
phi_angles, psi_angles = ramachandranangles(struc, stdresselector)
plot(x=map(rad2deg, phi_angles),
    y=map(rad2deg, psi_angles),
    Guide.xlabel("Phi / degrees"),
    Guide.ylabel("Psi / degrees"),
    Guide.xticks(ticks=[-180,-90,0,90,180]),
    Guide.yticks(ticks=[-180,-90,0,90,180]))
```

**F)** To calculate the RMSD and displacements between the heavy (non-hydrogen) atoms of two models in an NMR structure:

```julia
downloadpdb("1SSU")
struc_nmr = read("1SSU.pdb", PDB)
rmsd(struc_nmr[5], struc_nmr[10], heavyatomselector)
displacements(struc_nmr[5], struc_nmr[10], heavyatomselector)
```
