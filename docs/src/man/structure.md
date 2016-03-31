# Bio.Structure: Macromolecular Structures

    {meta}
    CurrentModule = Bio.Structure
    DocTestSetup = quote
        using Bio.Structure
        path = Pkg.dir("Bio", "test", "BioFmtSpecimens")
        if !isdir(path)
            run(`git clone --depth 1 https://github.com/BioJulia/BioFmtSpecimens.git $(path)`)
        end
        filepath_1EN2 = Pkg.dir("Bio", "test", "BioFmtSpecimens", "PDB", "1EN2.pdb")
    end

The `Bio.Structure` module provides functionality to read Protein Data Bank (PDB) files and manipulate macromolecular structures.


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

| Command                     | Returns                                                                         |
| :-------------------------- | :------------------------------------------------------------------------------ |
| `struc[1]`                  | Model 1                                                                         |
| `struc[1]['A']`             | Model 1, chain A                                                                |
| `struc['A']`                | The lowest model (model 1), chain A                                             |
| `struc['A']["50"]`          | Model 1, chain A, residue 50                                                    |
| `struc['A'][50]`            | Shortcut to above if it is not a hetero residue and the insertion code is blank |
| `struc['A']["H_90"]`        | Model 1, chain A, hetero residue 90                                             |
| `struc['A'][50]["CA"]`      | Model 1, chain A, residue 50, atom name CA                                      |
| `struc['A'][15]["CG"]['A']` | For disordered atoms, access a specific location                                |
| `struc['A'][15]["CG"]`      | For disordered atoms, access the default location                               |

Disordered atoms are stored in a `DisorderedAtom` container but calls fall back to the default atom, so disorder can be ignored if you are not interested in it.

Disordered residues (i.e. point mutations with different residue names) are stored in a `DisorderedResidue` container.

The idea is that disorder will only bother you if you want it to. See the [Biopython discussion](http://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ#How_is_disorder_handled.3F) for more.

Properties can be retrieved as follows:

| Command                  | Returns                                          |
| :----------------------- | :----------------------------------------------- |
| `structurename(struc)`   | Name of a structure                              |
| `modelnumbers(struc)`    | Model numbers in a structure                     |
| `chainids(struc)`        | Chain IDs in a structure                         |
| `modelnumber(model)`     | Number of a model                                |
| `chainids(model)`        | Chain IDs in a model                             |
| `chainid(chain)`         | Chain ID of a chain                              |
| `resids(chain)`          | Residue IDs in a chain                           |
| `resname(res)`           | Name of a residue                                |
| `chainid(res)`           | Chain ID of a residue                            |
| `resnumber(res)`         | Residue number of a residue                      |
| `inscode(res)`           | Insertion code of a residue                      |
| `ishetres(res)`          | `true` if the residue consists of hetero atoms   |
| `atomnames(res)`         | Atom names in a residue                          |
| `isdisorderedres(res)`   | `true` if the residue has multiple residue names |
| `resid(res)`             | Residue ID of a residue                          |
| `ishetatom(atom)`        | `true` if the atom is a hetero atom              |
| `serial(atom)`           | Serial number of an atom                         |
| `atomname(atom)`         | Name of an atom                                  |
| `altlocid(atom)`         | Alternative location ID of an atom               |
| `resname(atom)`          | Residue name of an atom                          |
| `chainid(atom)`          | Chain ID of an atom                              |
| `resnumber(atom)`        | Residue number of an atom                        |
| `inscode(atom)`          | Insertion code of an atom                        |
| `x(atom)`                | x coordinate of an atom                          |
| `y(atom)`                | y coordinate of an atom                          |
| `z(atom)`                | z coordinate of an atom                          |
| `coords(atom)`           | coordinates of an atom as a `Vector{Float64}`    |
| `occupancy(atom)`        | Occupancy of an atom (default is `1.0`)          |
| `tempfac(atom)`          | Temperature factor of an atom (default is `0.0`) |
| `element(atom)`          | Element of an atom (default is `""`)             |
| `charge(atom)`           | Charge of an atom (default is `""`)              |
| `isdisorderedatom(atom)` | `true` if the atom is disordered                 |
| `resid(atom)`            | Residue ID of an atom                            |

The coordinates of an atom can be set using `x!`, `y!`, `z!` and `coords!`.


## Manipulating Structures

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

`collect`, `collectresidues` and `collectatoms` can be used to get lists of sub-elements.

Selectors are functions passed as additional arguments to `collectresidues` and `collectatoms`. Only residues/atoms that return `true` when passed to the selector are retained.

| Command                                                 | Action                                                            |
| :------------------------------------------------------ | :---------------------------------------------------------------- |
| `collect(struc['A'][50])`                               | Collect the sub-elements of an element, e.g. atoms from a residue |
| `collectresidues(struc)`                                | Collect the residues of an element                                |
| `collectatoms(struc) `                                  | Collect the atoms of an element                                   |
| `collectatoms(struc, calphaselector)`                   | Collect the C-alpha atoms of an element                           |
| `collectatoms(struc, calphaselector, disorderselector)` | Collect the disordered C-alpha atoms of an element                |

It is easy to define your own selector. The below will collect all atoms with x coordinate less than 0:

```julia
xselector(atom::AbstractAtom) = x(atom) < 0
collectatoms(struc, xselector)
```

`countmodels`, `countchains`, `countresidues` and `countatoms` can be used to count elements. For example:

    {meta}
    DocTestSetup = quote
        using Bio.Structure
        filepath_1EN2 = Pkg.dir("Bio", "test", "BioFmtSpecimens", "PDB", "1EN2.pdb")
        struc = read(filepath_1EN2, PDB)
    end

```julia
julia> countatoms(struc)
754

julia> countatoms(struc, calphaselector)
85

julia> countresidues(struc, stdresselector)
85
```

`organise`, `organisemodel` and `organisestruc` can be used to organise sub-elements into elements:

| Command                                  | Action                                    |
| :--------------------------------------- | :---------------------------------------- |
| `organise(collectatoms(struc))`          | Organise an atom list into a residue list |
| `organise(collectresidues(struc))`       | Organise a residue list into a chain list |
| `organise(struc['A'])`                   | Organise chain(s) into a model            |
| `organise(struc[1])`                     | Organise model(s) into a structure        |
| `organisemodel(collectatoms(struc))`     | Organise elements into a model            |
| `organisestructure(collectatoms(struc))` | Organise elements into a structure        |

Distances can be calculated. The minimum distance between residue 10 and 20 is:

```julia
julia> distance(struc['A'][10], struc['A'][20])
10.782158874733762
```

RMSDs/displacements between elements of the same size can also be calculated with `rmsd` and `displacements`.


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
atoms = collectatoms(struc, calphaselector)
res_numbers = map(resnumber, atoms)
temp_facs = map(tempfac, atoms)
plot(x=res_numbers, y=temp_facs, Geom.line)
```

**B)** To find all C-alpha atoms within 5 Angstroms of residue 38:

```julia
for atom in atoms
    if distance(struc['A'][38], atom) < 5.0 && resnumber(atom) != 38
        show(atom)
    end
end
```

**C)** To calculate the RMSD and displacements between the heavy (non-hydrogen) atoms of two models in an NMR structure:

```julia
downloadpdb("1SSU")
struc_nmr = read("1SSU.pdb", PDB)
rmsd(struc_nmr[5], struc_nmr[10], heavyatomselector)
displacements(struc_nmr[5], struc_nmr[10], heavyatomselector)
```
