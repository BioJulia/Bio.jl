
# Parsing Phylogenetic trees from Newick formatted files.

module NewickParserImpl

import Bio.Ragel
import Bio
using BufferedStreams, Switch

immutable Newick <: Bio.FileFormat

%% {
    machine newick;

    action anchor { Ragel.@anchor! }
    action token {
        Ragel.@copy_from_anchor!(token)
        println(token)
    }

    length = '-'? digit+ ('.' digit*)?;
    unquoted_label = [A-Za-z_\-] [A-Za-z0-9_\-]*;
    quoted_label = '"' [^"]* '"';
    clade_entry = '(';
    clade_closure = ')';
    token = (clade_entry | clade_closure | ',' | ':' | ';' | length )

    main := space* (token space*)**;
} %%

%% write data;

type NewickParser <: Ragel.AbstractParser
    state::Ragel.State

    function NewickParser(input::BufferedInputStream)
        %% write init;

        return new(Ragel.State(cs, input))
    end
end

function Base.open(stream::BufferedInputStream, ::Type{Newick})
    return NewickParser(stream)
end

Ragel.@generate_read_function("newick", NewickParser, NewickTree,
    begin
        token = Bio.StringField()

        %% write exec;
    end)



end
