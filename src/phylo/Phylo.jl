#==================================#
#  Definition of the Phylo module  #
#==================================#

# Ben J. Ward, 2014.

module Phylo

using Docile

using Base.Intrinsics

using LightXML: XMLElement, get_elements_by_tagname, attribute

using Bio.Tools: Tokenizer, tokenize

import Base: getindex, setindex!, length, start, next, done, isempty, isequal, parent, delete!, search

import DataStructures: Queue, enqueue!, dequeue!, Stack, Queue, Deque

import Bio.FileFormat

export

  # PhyNodes and associated methods
  PhyNode, blisknown, confisknown, confidence, confidence!, isempty,
  name, branchlength, isleaf, haschildren, haschild, parentisself, hasparent,
  children, siblings, parent, isroot, isunlinked, islinked, isinternal,
  ispreterminal, issemipreterminal, countchildren, descendents,
  terminaldescendents, isancestral, mrca, name!, branchlength!,
  graft!, prune!, prunegraft!, delete!, detach!, isequal, distanceof,
  NodeCache, cachenodes,

  # Phylogeny and associated methods
  Phylogeny, isrooted, isrerootable, root!,
  setrerootable!, graft!, prune!, pruneregraft!, search,
  searchall, generateindex, ladderize, ladderize!,

  # PhylogenyIterator and associated methods
  PhylogenyIterator, DepthFirst, BreadthFirst, Tip2Root,
  getmrca, hasextensions, getroot, pathbetween,

  # File and format IO
  @newick_str, parsenewick, readnewick

include("nodes.jl")
include("phylogeny.jl")
include("treeio.jl")
include("iteration.jl")
include("annotation.jl")

end # module Phylo
