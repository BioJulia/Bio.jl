
using Bio.Intervals

#=

LiftOverChain is an implementation of a liftOver capable IntervalCollection

=#

immutable ChainBlock
   score::Float64
   tsize::Int64
   qname::String
   qsize::Int64
   qstart::Int64
   qend::Int64
   qstrand::Bio.Intervals.Strand
   id::Int64
   blocks::IntervalMap{Int64,Int64}
end

typealias LiftOverChain Bio.Intervals.IntervalCollection{ChainBlock}

# Replace this default io constructor with ragel parser?
function LiftOverChain(io)
   score     = 0.0
   tname     = ""
   tsize     = 0
   tstart    = 0
   tend      = 0
   tstrand   = '?'
   qname     = ""
   qsize     = 0
   qstart    = 0
   qend      = 0
   qstrand   = '?'
   id        = 0
   blocks  = IntervalMap{Int64,Int64}()

   chain = LiftOverChain()

   tcur   = 0
   qcur   = 0

   for l in eachline(io)
      spl = split(chomp(l), [' ','\t'])
      if l == "\n"
         # ignore
      elseif length(spl) == 13 && spl[1] == "chain"
         # header
         score = parse(Float64, spl[2])
         tname,tsize,tstrand = spl[3], parse(Int, spl[4]), spl[5][1] 
         tstart,tend         = parse(Int, spl[6]), parse(Int, spl[7] )
         qname,qsize,qstrand = spl[8], parse(Int, spl[9]), spl[10][1]
         qstart,qend         = parse(Int, spl[11]), parse(Int, spl[12])
         id                  = parse(Int, spl[13])
         blocks = IntervalTrees.IntervalMap{Int64,Int64}()
         tcur,qcur = tstart,qstart
      elseif length(spl) == 3
         # chain block
         blocksize  = parse(Int, spl[1])
         tinc, qinc = parse(Int, spl[2]), parse(Int, spl[3])          
         blocks[(tcur,tcur+blocksize)] = qcur
         tcur += blocksize + tinc
         qcur += blocksize + qinc
      elseif length(spl) == 1
         # last block
         blocksize = parse(Int, spl[1])
         blocks[(tcur,tcur+blocksize)] = qcur
         block = ChainBlock(score, tsize, String(qname), qsize, qstart, qend, 
                            Bio.Intervals.Strand(Char(qstrand)), id, blocks)
         targ  = Bio.Intervals.Interval(String(tname), tstart, tend, Char(tstrand), block)
         push!(chain, targ)
      else
         error("Malformed chain input file!")
      end
   end
   chain
end

function liftover{T}( chain::LiftOverChain, istream::Bio.Intervals.IntervalStreamOrArray{T}; minidentity=0.95 )
   lifted  = Vector{Nullable{Bio.Intervals.Interval{T}}}()

   cname   = ""
   cfirst  = 0
   clast   = 0
   cstrand = Bio.Intervals.STRAND_NA

   overlap = 0

   chain_state   = start(chain)
   istream_state = start(istream)

   if !done(chain, chain_state) && !done(istream, istream_state)
      chain_el   = next(chain, chain_state)
      istream_el = next(istream, istream_state)
   end

   while !done(chain, chain_state)

      if chain_state == start(chain) && istream_state == start(istream)
         chain_el,   chain_state   = next(chain, chain_state)
         istream_el, istream_state = next(istream, istream_state)
      end

      if precedes( chain_el, istream_el, Bio.Intervals.alphanum_isless ) 
         done(chain, chain_state) && break
         # increment chain
         chain_el, chain_state = next(chain, chain_state)
         cname   = chain_el.metadata.qname
         cstrand = chain_el.metadata.qstrand

      elseif precedes( istream_el, chain_el, Bio.Intervals.alphanum_isless )
         # push nullable{}
         push!( lifted, Nullable{Bio.Intervals.Interval{T}}() )
         done(istream, istream_state) && break
         # increment istream
         istream_el, istream_state = next(istream, istream_state)

      else # chain entry and interval overlap, now liftover

         block       = chain_el.metadata.blocks
         block_state = start(block)
         if !done(block, block_state)
            block_el = next(block, block_state)
         end

         while !done(block, block_state) 

            if block_state == start(block)
               block_el, block_state = next(block, block_state)
            end   

            if block_el.last < istream_el.first
               # increment block
               block_el, block_state = next(block, block_state)

            elseif istream_el.last < block_el.first
               # push empty 'deleted' coordinate as Nullable{I}
               push!( lifted, Nullable{Bio.Intervals.Interval{T}}() )
               # increment istream
               istream_el, istream_state = next(istream, istream_state)

            else
               # there is overlap with the current block
               overlap = istream_el.last - istream_el.first # initial

               if istream_el.first < block_el.first
                  # first coordinate is deleted.
                  cfirst   = block_el.value
                  overlap -= block_el.first - istream_el.first
               else
                  cfirst   = liftinternal( block_el, istream_el.first )
               end
               clast = liftsecond( block, block_el, block_state, istream_el )
               interval = Bio.Intervals.Interval(cname, cfirst, clast, cstrand, istream_el.metadata)
               push!( lifted, Nullable(interval) )
               done(istream, istream_state) && return lifted
               istream_el, istream_state = next(istream, istream_state)
               break
            end

         end # end while 
      end

   end # end while
   lifted
end

# private function called from `liftover` to lift the end coordinate
# of an interval with a copied version of the block iterator
function liftsecond( block, orig_el, orig_state, istream_el )

   block_state = deepcopy(orig_state)
   block_el    = orig_el
   prev_el     = nothing

   while !done(block, block_state)

      if prev_el != nothing
         block_el, block_state = next(block, block_state)
      end

      if block_el.first <= istream_el.last <= block_el.last
         # liftover coordinate
         return liftinternal( block_el, istream_el.last )
      elseif istream_el.last < block_el.first
         # last coordinate is deleted
         if prev_el == orig_el
            error("End coordinate cannot be less than the first!")
         else
            return liftlast( prev_el )
         end
      end
      prev_el = block_el
   end
   # hit the end of the chain without lifting
   # so return the end of the chain otherwise
   # we would need to handle multi-coordinate
   # liftover that could lift one interval to
   # multiple chrom-scaffold/strand
   return liftlast( block_el )
end

liftfirst( i::IntervalValue{Int64,Int64} ) = i.value
liftlast( i::IntervalValue{Int64,Int64} ) = i.value + (i.last - i.first)

function liftinternal( i::IntervalValue{Int64,Int64}, tolift::Int64 )
   if i.first <= tolift <= i.last
      return i.value + (tolift - i.first)
   else
      error( "Cannot lift coordinate that is not within block!" )
   end
end

function precedes{S,T}(a::Bio.Intervals.Interval{S}, b::Bio.Intervals.Interval{T}, func::Function)
    return (a.last < b.first && a.seqname == b.seqname) ||
        Bio.Intervals.alphanum_isless(a.seqname, b.seqname)::Bool
end

function Base.collect( chain::LiftOverChain )
   # this function could return the LiftOverChain as 
   # a vector of Interval{Interval}() pairs
   # where each pair describes a chain block in the
   # form of `from->to`
end
