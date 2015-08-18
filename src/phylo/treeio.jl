#==========================================#
# IO phylogenies stored in various formats #
#==========================================#

# Exception types
# ----------------

"""
A simple Exception type that is thrown by newick related functions when an error occurs.

**Fields:**

* `msg`: A `String` containing the message to print to screen with `showerror`.
"""
type NewickException <: Exception
    msg::ASCIIString
end

"""
Basic function that prints NewickExceptions to screen.
"""
Base.showerror(io::IO, e::NewickException) = print(io, "Error parsing newick string: ", e.msg)

# Internal methods used by newick string parser method.

"""
Makes a new clade when building a tree from a newick string.
This method of the function accepts no parameters and returns an empty `PhyNode`.
"""
function makenewclade()
  return PhyNode()
end

"""
Makes a new clade when building a tree from a newick string.

This method is used in the `parsenewick` method to take car of linking a
newly created node to its parent on creation.

**Parameters:**

* `parent`: The parent of the to-be-created PhyNode.

**Returns:** A reference to the newly created `PhyNode`.
"""
function makenewclade(parent::PhyNode)
    newclade = PhyNode()
    graft!(parent, newclade)
    return newclade
end

"""
Finishes the processing of the current clade in the newick file.

**Parameters:**

* `node`:            The `PhyNode` to finish processing.
* `valuesareconf`:   `Bool` that specifies whether the values of the clade are confidence values.
* `commentsareconf`: `Bool` that specifies if comments are confidence values.

**Returns:** The parent `PhyNode` of the `PhyNode` provided as the `node` parameter.
"""
function processclade(node::PhyNode, valuesareconf::Bool, commentsareconf::Bool)
    # Check if the node has a name, and if values are not confidence, and there are no conf
    # values in values or comments, and confience is not known.
    if name(node) != "" && !(valuesareconf || commentsareconf) && !confisknown(node) && haschildren(node)
        confidence!(node, parseconfidence(name(node)))
        if confisknown(node)
            name!(node, "")
        end
    end
    if hasparent(node)
        parent = node.parent
        prune!(node)
        graft!(parent, node)
        return parent
    end
end

"""
Parses confidence values from string - if this is not possible,
-1.0 is returned, which is the value for 'unknown' in `Phylo`.
"""
function parseconfidence(text::ASCIIString)
    try
        return float(text)
    catch
        return -1.0
    end
end

# Method for parsing a Newick formatted string.

"""
Build a `Phylogeny` from a String formatted as a Newick string.

Newick strings are of the form `(((A,B),C),D);`.
In such strings, parentheses delimit clades, and text delimits taxa names. Somtimes accompanied
by a floating point value that may be a branch length, or clade support value.

**Paramerters:**

* `newickstring`: A newick formatted `String`.
* `commentsareconf`: `Bool` value indicating if comments provide clade support values.
* `valuessareconf`: `Bool` value indicating if values (usually branchlengths) provide clade support values.

**Returns:** A `Phylogeny` constructed from the string.
"""
function parsenewick(newickstring::ASCIIString, commentsareconf::Bool = false, valuesareconf::Bool = false)
  # Create a definition of the tokens that appear in a newick string, and the meanings of them.
    definition::Vector{Tuple{ASCIIString, Regex}} = [("open paren", r"\("),
                                                    ("close paren", r"\)"),
                                                    ("unquoted node label", r"[^\s\(\)\[\]\'\:\;\,]+"),
                                                    ("edge length", r"\:[0-9]*\.?[0-9]+([eE][+-]?[0-9]+)?"),
                                                    ("comma", r"\,"),
                                                    ("comment", r"\[(\\.|[^\]])*\]"),
                                                    ("quoted node label", r"\'(\\.|[^\'])*\'"),
                                                    ("semicolon", r"\;"), ("newline", r"\n")]
    tokenizer::Tokenizer = Tokenizer(definition)
    # Convet the newick string into a series of tokens than can be considered in turn and understood.
    tokens = tokenize(strip(newickstring), tokenizer).tokens
    # Create the first clade, i.e. the root and set the variable that points to the current clade
    # to the root.
    root = PhyNode("Root")
    phy = Phylogeny("", root, true, true)
    current = root
    @assert root === current
    enteringbl = false
    # Assign two variables which will keep track of the number of left and right parentheses.
    leftpcount = 0
    rightpcount = 0
    state = start(tokens)
    comments = TreeAnnotations(phy, String)
    while !done(tokens, state)
        token, state = next(tokens, state)
        # TODO - Do something about this truly horrible set of chained if statements.
        if startswith(token, "\'")
            # This is a quoted label, characters need to be added to the clade name.
            name!(current, name(current) * token[2:end])
        elseif startswith(token, "[")
            # This is a comment for the clade.
            comments[current] = token[2:end-1]
            if commentsareconf
                confidence!(current, parseconfidence(comments[current]))
            end
        elseif token == "("
            # The start of a new clade. It is a child of the current clade.
            current = makenewclade(current)
            enteringbl = false
            leftpcount += 1
        elseif token == ","
            # If the current clade is the root, it means the external parentheses are missing.
            if current === root
                root = makenewclade()
                name!(root, "Root")
                prune!(current)
                graft!(root, current)
            end
            # Start a new child clade at the same level as the current clade. i.e. a sibling.
            parent = processclade(current, valuesareconf, commentsareconf)
            current = makenewclade(parent)
            enteringbl = false
        elseif token == ")"
            # Addition of children for the current clade is done, so we process it.
            parent = processclade(current, valuesareconf, commentsareconf)
            parent === current ? throw(NewickException("Tried to ascend to a parent and failed.")) : nothing
            current = parent
            enteringbl = false
            rightpcount += 1
        elseif token == ";"
            break
        elseif startswith(token, ":")
            # This token is branchlength or confidence...
            value = float(token[2:end])
            if valuesareconf
                confidence!(current, value)
            else
                branchlength!(current, value)
            end
        elseif token == "\n"
            continue
        else
            # Current token is an unquoted node label, and so we simply need to set the clade
            # name to the token.
            name!(current, token)
        end
    end
    if leftpcount != rightpcount
        throw(NewickException("The number of left and right parentheses do not match."))
    end
    try
        token, state = next(tokens, state)
        throw(NewickException("There are characters after the semicolon in the newick string."))
    catch exception
        if isa(exception, NewickException)
            rethrow(exception)
        end
    end
    processclade(current, valuesareconf, commentsareconf)
    processclade(root, valuesareconf, commentsareconf)
    return phy, comments
end

# Method for reading a file of Newick formatted strings i.e. Newick Format files.

function readnewick(file::ASCIIString)
    instream = open(expanduser(file))
    instring = readall(instream)
    close(instream)
    newickstrings = split(instring, ';')
    newickstrings = [replace(i, r"(\r|\n|\s)", "") for i in newickstrings]
    newickstrings = newickstrings[bool([length(t) > 0 for t in newickstrings])]
    if length(newickstrings) > 1
        return [parsenewick(i) for i in newickstrings]
    else
        return parsenewick(newickstrings[1])
    end
    error("No phylogenies detected in the file.")
end





# Recusrsiveley builds tree structure from PhyloXML format.
# function buildrecursively(xmlclade::XMLElement, extensionsArray)
#   label::ASCIIString = ""
#   node::PhyNode = PhyNode(name = label, branchlength = (bl == nothing ? -1.0 : float(bl)), ext = PhyExtension[parsephylo(xmlclade, i) for i in extensionsArray])
#   xmlchildren::Array{XMLElement, 1} = get_elements_by_tagname(xmlclade, "clade")
#   for n in xmlchildren
#     graft!(node, buildrecursively(n, extensionsArray))
#   end
#   return node
# end

# function parsephyloxml(xmltree::XMLElement, ::Type{Phylogeny}, extensions) # extensions is a tuple of types.
#   # Get the first clade of the tree - the root clade in the XML, the recursiveley build up the tree structure.
#   # Function assumes a tree is rerootable unless otherwise specified, and assumes it's not rooted unless otherwise specified in the file.
#   # Function also assumes that there is and can only be one name for the tree.
#   XML::XMLElement = get_elements_by_tagname(xmltree, "clade")[1]
#   root::PhyNode = buildRecursively(XML, extensions, nothing)
#   name::Array{XMLElement, 1} = get_elements_by_tagname(xmltree, "name")
#   treename = length(name) == 1 ? content(name[1]) : ""
#   rerootable::Bool = attribute(xmltree, "rerootable") == "false" ? true : false
#   rooted::Bool = attribute(xmltree, "rooted") == "true" ? true : false
#   return Phylogeny(treename, root, rooted, rerootable)
# end

# # Open to suggestions making this function neater re deciding on different formats. Since branching conditionally on the format variable seems less neat to me,
# # than taking advantage of method dispatch. Perhaps some dummy types for the input formats to allow method dispatch would be better?

# function readphyloxml(file::ASCIIString, Extensions...)
#   treedoc::XMLDocument = parse_file(expanduser(file))
#   phylogenies::Array{XMLElement,1} = get_elements_by_tagname(root(treedoc), "phylogeny")
#   if length(phylogenies) > 1
#     return [parsephyloxml(i, Phylogeny, Extensions) for i in phylogenies]
#   else
#     return parsephyloxml(phylogenies[1], Phylogeny, Extensions)
#   end
#   error("No Trees found in file...")
# end


### Newick writing functionality...

function informationmaker(plain::Bool, blconfidence::Bool, blonly::Bool,
    maximumconf::Float64)
    if plain
        #If plain is selected, a plain newick string is to be made
        makeinfostring(node, terminal = false) = "$(comment(node))"
    elseif blconfidence
        # Branchlengths are support values. Therefore you want to ignore
        # the actual branchlengths, instead writing out the confidence to the string.
        function makeinfostring(node, terminal = false)
            if terminal
                # Terminal branches by definition, have 100% support as they are what
                # is observed in nature.
                out = ":$(maximumconf)$(comment(node))"
            else
                out = ":$(confidence(node))$(comment(node))"
            end
            return out
        end
    elseif blonly
        # The function should ignore confidence and only write branchlength.
        function makeinfostring(node, terminal = false)
            return ":$(branchlength(node))$(comment(node))"
        end
    else
        function makeinfostring(node, terminal = false)
            if terminal || !hasconfidence(node)
                return ":$(branchlength(node))$(comment(node))"
            else
                return "$(confidence(node):$(branchlength(node))$(comment(node)))"
            end
        end
    end
    return makeinfostring
end

function newickstring(blconfidence = false, blonly = false, plain = false,
    plainnewick = true, ladderize = nothing, maximumconf = 1.0, confSF = "1.2f",
    formatbl = "1.5f")

    if blconfidence || blonly
        plain = false
    end

    infostring = informationmaker(plain, blconfidence, blonly,
        maximumconf, formatconf, formatbl)
end
