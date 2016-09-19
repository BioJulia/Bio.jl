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

The `Bio.Structure` module provides functionality to manipulate macromolecular structures, and in particular to read and write [Protein Data Bank](http://www.rcsb.org/pdb/home/home.do) (PDB) files. It is designed to be used for standard structural analysis tasks, as well as acting as a platform on which others can build to create more specific tools.


## Parsing PDB files

To download a PDB file:

```julia
downloadpdb("1EN2")
```

To parse a PDB file into a Structure-Model-Chain-Residue-Atom framework:

```julia
julia> struc = read(filepath_1EN2, PDB)
Bio.Structure.ProteinStructure
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

| Command                     | Returns                                                                         | Return type       |
| :-------------------------- | :------------------------------------------------------------------------------ | :---------------- |
| `struc[1]`                  | Model 1                                                                         | `Model`           |
| `struc[1]['A']`             | Model 1, chain A                                                                | `Chain`           |
| `struc['A']`                | The lowest model (model 1), chain A                                             | `Chain`           |
| `struc['A']["50"]`          | Model 1, chain A, residue 50                                                    | `AbstractResidue` |
| `struc['A'][50]`            | Shortcut to above if it is not a hetero residue and the insertion code is blank | `AbstractResidue` |
| `struc['A']["H_90"]`        | Model 1, chain A, hetero residue 90                                             | `AbstractResidue` |
| `struc['A'][50]["CA"]`      | Model 1, chain A, residue 50, atom name CA                                      | `AbstractAtom`    |
| `struc['A'][15]["CG"]['A']` | For disordered atoms, access a specific location                                | `Atom`            |

Disordered atoms are stored in a `DisorderedAtom` container but calls fall back to the default atom, so disorder can be ignored if you are not interested in it.

Disordered residues (i.e. point mutations with different residue names) are stored in a `DisorderedResidue` container.

The idea is that disorder will only bother you if you want it to. See the [Biopython discussion](http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ#How_is_disorder_handled.3F) for more.

Properties can be retrieved as follows:

| Function           | Returns                                                       | Return type                     |
| :----------------- | :------------------------------------------------------------ | :------------------------------ |
| `serial`           | Serial number of an atom                                      | `Int`                           |
| `atomname`         | Name of an atom                                               | `String`                        |
| `altlocid`         | Alternative location ID of an atom                            | `Char`                          |
| `x`                | x coordinate of an atom                                       | `Float64`                       |
| `y`                | y coordinate of an atom                                       | `Float64`                       |
| `z`                | z coordinate of an atom                                       | `Float64`                       |
| `coords`           | coordinates of an atom                                        | `Array{Float64,1}`              |
| `occupancy`        | Occupancy of an atom (default is `1.0`)                       | `Float64`                       |
| `tempfactor`       | Temperature factor of an atom (default is `0.0`)              | `Float64`                       |
| `element`          | Element of an atom (default is `"  "`)                        | `String`                        |
| `charge`           | Charge of an atom (default is `"  "`)                         | `String`                        |
| `residue`          | Residue an atom belongs to                                    | `Residue`                       |
| `ishetero`         | `true` if the residue or atom is a hetero residue/atom        | `Bool`                          |
| `isdisorderedatom` | `true` if the atom is disordered                              | `Bool`                          |
| `resname`          | Residue name of a residue or atom                             | `String`                        |
| `resnumber`        | Residue number of a residue or atom                           | `Int`                           |
| `inscode`          | Insertion code of a residue or atom                           | `Char`                          |
| `resid`            | Residue ID of an atom or residue (`full=true` includes chain) | `String`                        |
| `atomnames`        | Atom names of the atoms in a residue, sorted by serial        | `Array{String,1}`               |
| `atoms`            | Dictionary of atoms in a residue                              | `Dict{String, AbstractAtom}`    |
| `isdisorderedres`  | `true` if the residue has multiple residue names              | `Bool`                          |
| `disorderedres`    | Access a particular residue name in a `DisorderedResidue`     | `Residue`                       |
| `chain`            | Chain a residue or atom belongs to                            | `Chain`                         |
| `chainid`          | Chain ID of a chain, residue or atom                          | `Char`                          |
| `resids`           | Sorted residue IDs in a chain                                 | `Array{String,1}`               |
| `residues`         | Dictionary of residues in a chain                             | `Dict{String, AbstractResidue}` |
| `model`            | Model a chain, residue or atom belongs to                     | `Model`                         |
| `modelnumber`      | Model number of a model, chain, residue or atom               | `Int`                           |
| `chainids`         | Sorted chain IDs in a model or structure                      | `Array{Char,1}`                 |
| `chains`           | Dictionary of chains in a model or structure                  | `Dict{Char, Chain}`             |
| `structure`        | Structure a model, chain, residue or atom belongs to          | `ProteinStructure`              |
| `structurename`    | Name of the structure an element belongs to                   | `String`                        |
| `modelnumbers`     | Sorted model numbers in a structure                           | `Array{Int,1}`                  |
| `models`           | Dictionary of models in a structure                           | `Dict{Int, Model}`              |

The `spaces` keyword argument determines whether surrounding whitespace is retained for `atomname`, `element`, `charge`, `resname` and `atomnames` (default `false`).

The coordinates of an atom can be set using `x!`, `y!`, `z!` and `coords!`.


## Manipulating structures

Elements can be looped over to reveal the sub-elements in the correct order:

```julia
for mod in struc
    for ch in mod
        for res in ch
            for at in res
                # Do something
            end
        end
    end
end
```

Models are ordered numerically; chains are ordered by character, except the empty chain is last; residues are ordered by residue number and insertion code with hetero residues after standard residues; atoms are ordered by atom serial.

`collect` can be used to get arrays of sub-elements. `collectatoms`, `collectresidues`, `collectchains` and `collectmodels` return arrays of a particular type from a structural element or element array.

Selectors are functions passed as additional arguments to these functions. Only elements that return `true` when passed to the selector are retained.

| Command                                                 | Action                                                            | Return type                |
| :------------------------------------------------------ | :---------------------------------------------------------------- | :------------------------- |
| `collect(struc['A'][50])`                               | Collect the sub-elements of an element, e.g. atoms from a residue | `Array{AbstractAtom,1}`    |
| `collectresidues(struc)`                                | Collect the residues of an element                                | `Array{AbstractResidue,1}` |
| `collectatoms(struc) `                                  | Collect the atoms of an element                                   | `Array{AbstractAtom,1}`    |
| `collectatoms(struc, calphaselector)`                   | Collect the C-alpha atoms of an element                           | `Array{AbstractAtom,1}`    |
| `collectatoms(struc, calphaselector, disorderselector)` | Collect the disordered C-alpha atoms of an element                | `Array{AbstractAtom,1}`    |

It is easy to define your own atom, residue, chain or model selectors. The below will collect all atoms with x coordinate less than 0:

```julia
xselector(at::AbstractAtom) = x(at) < 0
collectatoms(struc, xselector)
```

Alternatively, you can use an anonymous function:

```julia
collectatoms(struc, at -> x(at) < 0)
```

`countatoms`, `countresidues`, `countchains` and `countmodels` can be used to count elements. For example:

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

julia> countresidues(struc, standardselector)
85
```

The sequence of a protein can be retrieved by passing a `Chain` or array of residues to `AminoAcidSequence`:

```julia
julia> AminoAcidSequence(struc['A'], standardselector)
85aa Amino Acid Sequence:
RCGSQGGGSTCPGLRCCSIWGWCGDSEPYCGRTCENKCWSGERSDHRCGAAVGNPPCGQDRCCSVHGWCGGGNDYCSGGNCQYRC
```


## Writing PDB files

PDB format files can be written:

```julia
writepdb("1EN2_out.pdb", struc)
```

Any element type can be given as input to `writepdb`. Atom selectors can also be given as additional arguments:

```julia
writepdb("1EN2_out.pdb", struc, backboneselector)
```


## Spatial calculations

Various functions are provided to calculate spatial quantities on proteins:

| Command              | Returns                                                                                         |
| :------------------- | :---------------------------------------------------------------------------------------------- |
| `distance`           | Minimum distance between two elements                                                           |
| `sqdistance`         | Minimum square distance between two elements                                                    |
| `bondangle`          | Angle between three atoms                                                                       |
| `dihedralangle`      | Dihedral angle defined by four atoms                                                            |
| `omegaangle`         | Omega angle between residue and previous residue                                                |
| `phiangle`           | Phi angle between residue and previous residue                                                  |
| `psiangle`           | Psi angle between residue and next residue                                                      |
| `ramachandranangles` | `Vector`s of phi and psi angles of an element                                                   |
| `contactmap`         | Contact map of two element, or one element with itself                                          |
| `rmsd`               | RMSD between two elements of the same size - assumes they are superimposed                      |
| `displacements`      | `Vector` of displacements between two elements of the same size - assumes they are superimposed |

For example:

```julia
julia> distance(struc['A'][10], struc['A'][20])
10.782158874733762

julia> rad2deg(bondangle(struc['A'][50]["N"], struc['A'][50]["CA"], struc['A'][50]["C"]))
110.77765846083398

julia> rad2deg(dihedralangle(struc['A'][50]["N"], struc['A'][50]["CA"], struc['A'][50]["C"], struc['A'][51]["N"]))
-177.38288114072924

julia> rad2deg(psiangle(struc['A'][50], struc['A'][51]))
-177.38288114072924
```


## Examples

A few further examples of `Bio.Structure` usage are given below.

**A)** To plot the temperature factors of a protein, if you have Gadfly installed:

```julia
using Gadfly
calphas = collectatoms(struc, calphaselector)
plot(x=resnumber.(calphas),
     y=tempfactor.(calphas),
     Guide.xlabel("Residue number"),
     Guide.ylabel("Temperature factor"),
     Geom.line)
```

**B)** To print the PDB records for all C-alpha atoms within 5 Angstroms of residue 38:

```julia
for at in calphas
    if distance(struc['A'][38], at) < 5.0 && resnumber(at) != 38
        println(join(pdbline(at)))
    end
end
```

**D)** To view the contact map of a structure:

```julia
cbetas = collectatoms(struc, cbetaselector)
contacts = contactmap(cbetas, 7.0)
for i in 1:length(cbetas)
    for j in 1:length(cbetas)
        if contacts[i,j]
            print("█")
        else
            print(" ")
        end
    end
    println()
end
```

`cbetaselector` selects C-beta atoms, or C-alpha atoms for glycine residues. `contactmap` can also be given two structural elements as arguments, in which case a non-symmetrical 2D array is returned showing contacts between the elements.

**E)** To show the Ramachandran phi/psi angle plot of a structure, if you have Gadfly installed:

```julia
using Gadfly
phi_angles, psi_angles = ramachandranangles(struc, standardselector)
plot(x=rad2deg.(phi_angles),
     y=rad2deg.(psi_angles),
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
