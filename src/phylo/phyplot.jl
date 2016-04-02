@enum nodesym circle=1 square=2 rect=3
@enum textpos auto=1 alignright=2 alignleft=3 above=4
@enum plottype phylogram=1 fan=2

type NodeFormat#{S<:AbstractString}
    #fill::NodeAnnotations{Colorant}
    #stroke::NodeAnnotations{Colorant}
    #strokewidth::NodeAnnotations{Measure}
    #label_above::NodeAnnotations{S}
    #label_below::NodeAnnotations{S}
    #symbol::nodesym
    #symbol_size::NodeAnnotations{Float64}
    #align_label_above::textpos
    #align_label_below::textpos
end

type EdgeFormat
    line_width::Float64
    line_color::Union{Colorant, NodeAnnotations{Colorant}}
end

function EdgeFormat(;linewidth::Float64 = 1., linecolor::Union{Colorant, NodeAnnotations{Colorant}} = colorant"black")
    EdgeFormat(linewidth, linecolor)
end

type TipFormat#{S<:AbstractString}
    font_size::Float64
    #labels::NodeAnnotations{S} #This is meant to be the names dict when that gets implemented
    #alignment::textpos
    #symbol::nodesym
    #symbol_stroke::NodeAnnotations{Colorant}
end

TipFormat(;fontsize = 1.) = TipFormat(fontsize)

function termdec(phy::Phylogeny)
    termdec(phy.root)
end

function termdec(node::PhyNode) #replacements for terminaldescendents which isn't working
    ret = Array{PhyNode, 1}(0)
    for i in DepthFirst(node)
        if isleaf(i)
            push!(ret, i)
        end
    end
    ret
end


function phyplot(phy::Phylogeny, plot_type::plottype=phylogram; edges::EdgeFormat = EdgeFormat(), nodes::NodeFormat = NodeFormat(), tips::TipFormat = TipFormat(), rootbranch::Bool = (plot_type == fan)) #replace with a kwargs sloop like in GadFly

    function findxy(phy::Phylogeny)
        height = [tip => float(i) for (i, tip) in enumerate(termdec(phy))]

        function loc(clade::PhyNode)
            if !in(clade, keys(height))
                for subclade in children(clade)
                    loc(subclade)
                end
            end
            if !isleaf(clade)
                ch_heights = [height[child] for child in children(clade)]
                height[clade] = (maximum(ch_heights) + minimum(ch_heights)) / 2.
            end
        end
        loc(phy.root)

        depth = NodeAnnotations{Float64}(phy.root => 0)
        for(clade in DepthFirst(phy))
            if !parentisself(clade)
                depth[clade] = depth[parent(clade)] + branchlength(clade, 1.)
            end
        end
        depth, height
    end

    function compose_segments(x1::Float64, y1::Float64, x2::Float64, y2::Float64)
        [(x1, y1), (x2, y2)]
    end

    function process_edges!(ed::EdgeFormat, x::Phylogeny)
        if isa(ed.line_color, Colorant)
            col::Colorant = ed.line_color
            ed.line_color = [node => col for node in DepthFirst(phy)]
        end
    end

    process_edges!(edges, phy)

    x, y = findxy(phy)

    if rootbranch
        for key in keys(x)
            x[key] += maximum(values(x)) * 0.05
        end
    end

    internals = setdiff(keys(x), termdec(phy))
    horizontal_lines = [compose_segments(x[node], y[node], parentisself(node) ? 0. : x[parent(node)], y[node]) for node in keys(x)]
    vertical_lines1 = [compose_segments(x[node], y[children(node)[1]], x[node], y[node]) for node in internals]
    vertical_lines2 = [compose_segments(x[node], y[node], x[node], y[children(node)[end]]) for node in internals]
    horiz_color = [edges.line_color[node] for node in keys(x)]
    vert1_color = [edges.line_color[children(node)[1]] for node in internals]
    vert2_color = [edges.line_color[children(node)[2]] for node in internals]

    phyplot_intern(Val{plot_type}, x, y, horizontal_lines, vertical_lines1, vertical_lines2, horiz_color, vert1_color, vert2_color, edges, tips, termdec(phy))
end

function phyplot_intern(::Type{Val{phylogram}}, x, y, horizontal_lines, vertical_lines1, vertical_lines2, horiz_color, vert1_color, vert2_color, edges, tips, terminals)
    xtext = Float64[]
    ytext = Float64[]
    namestext = UTF8String[]

    for sp in terminals
        push!(xtext, x[sp]*1.02)
        push!(ytext, y[sp] + 0.5)
        push!(namestext, name(sp))
    end

    # # code for showing only some species
    # specspace = maximum(values(y)) > show_tips ? div(Int(maximum(values(y))), show_tips) : 1
    # showspecs = specspace:specspace:Int(maximum(values(y)))
    # xtext = xtext[showspecs]
    # ytext = ytext[showspecs]
    # namestext = namestext[showspecs]

    compose(
    context(units = UnitBox(-0.02*maximum(values(x)),-0.02 * minimum(values(y)),maximum(values(x))*1.4, maximum(values(y))*1.1)),
            (context(), line(horizontal_lines), stroke(horiz_color)),
            (context(), line(vertical_lines1), stroke(vert1_color)),
            (context(), line(vertical_lines2), stroke(vert2_color)),
            (context(), text(xtext, ytext, namestext),
    fontsize(minimum([3.5, tips.font_size * 120/length(terminals)]))),  linewidth(edges.line_width * 0.2)
    )
end

function phyplot_intern(::Type{Val{fan}}, x, y, horizontal_lines, vertical_lines1, vertical_lines2, horiz_color, vert1_color, vert2_color, edges, tips, terminals)

    function anglefrompoint(point::Tuple{Number, Number})
        point[2] >= 0 ? acos(point[1]) : pi + acos(-point[1])
    end

    function beziercurve(start::Tuple{Number, Number}, angle::Number, center::Tuple{Number, Number})
        pt = start[1] - center[1], start[2] - center[2]
        ret = beziercurve(pt, angle)
        [(pt[1] + center[1], pt[2] + center[2]) for pt in ret]
    end

    function beziercurve(start::Tuple{Number, Number}, angle::Number)
        r = sqrt(start[1]^2 + start[2]^2)
        start = start[1]/r, start[2]/r
        arc = beziercurve(angle)
        ret = [rotatepoint(pt, anglefrompoint(start)) for pt in arc]
        [(r*pt[1], r*pt[2]) for pt in ret]
    end

    function beziercurve(angle::Number)
        L  = 4//3*tan(angle/4.)
        p1 = 1., 0.
        p2 = 1., L
        p3 = cos(angle) + L*sin(angle), sin(angle) - L*cos(angle)
        p4 = cos(angle), sin(angle)
        [p1, p2, p3, p4]
    end

    function rotatepoint(point::Tuple{Number, Number}, angle::Number)
        point[1]*cos(angle) - point[2]*sin(angle), point[1]*sin(angle) + point[2]*cos(angle)
    end

    function pointcoords(angle::Float64, radius::Float64 = 1.)
        radius*sin(angle), radius*cos(angle)
    end

    function normpoint(pt::Tuple{Float64, Float64})
        pointcoords(pt[2] * circval, pt[1] * radval)
    end

    function normpoint(ar::Array{Tuple{Float64, Float64},1})
        [normpoint(x) for x in ar]
    end

    circval = 2pi/(1 + maximum(values(y)))
    radval = 1 / maximum(values(x))

    depth = x
    height = y

    elements = Any[]

    for lin in vertical_lines1
        bez = beziercurve(normpoint(lin[1]), (lin[1][2] - lin[2][2]) * circval)
        push!(elements, compose(context(),curve(bez...), stroke(vert1_color)))
    end

    for lin in vertical_lines2
        bez = beziercurve(normpoint(lin[1]), (lin[1][2] - lin[2][2]) * circval)
        push!(elements, compose(context(),curve(bez...), stroke(vert2_color)))
    end

    rad = [normpoint(x) for x in horizontal_lines]
    push!(elements, compose(context(),line(rad), stroke(horiz_color)))

    haligns = HAlignment[]
    valigns = VAlignment[]
    txts = UTF8String[]
    rots = Rotation[]
    pointx = Float64[]
    pointy = Float64[]

    for tip in terminals
        point = normpoint((1.03*depth[tip], height[tip]))
        push!(pointx, point[1])
        push!(pointy, point[2])
        push!(valigns, VCenter())
        push!(haligns, point[1] < 0 ? HRight() : HLeft())
        push!(txts, tip.name)
        rot = point[1] < 0 ?  pi + pi/2. - height[tip]*circval : pi/2. - height[tip]*circval
        push!(rots, Rotation(rot, point))
    end

    #compose(context(0,0, min(h, w), min(h, w), units = UnitBox(-1.6, -1.6, 3.2, 3.2)),
    compose(context(0,0, h, h, units = UnitBox(-1.6, -1.6, 3.2, 3.2)),
        elements...,
        text(pointx, pointy, txts, haligns, valigns, rots),
    fontsize(minimum([3, tips.font_size * 200/length(terminals)])),  linewidth(edges.line_width * 0.2)
    )
end
