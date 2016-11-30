# EUtils
# ======
#
# APIs for E-Utilities.
#
# This file is a part of BioJulia.
# License is MIT: https://github.com/BioJulia/Bio.jl/blob/master/LICENSE.md

"""
Entrez Programming Utilities (or E-Utilities) module.

The `EUtils` module provides a programming interface to the Entrez databases at
NCBI. Nine functions are exported from this module:
1. `einfo`: database statistics
2. `esearch`: text search
3. `epost`: UID upload
4. `esummary`: document summary download
5. `efetch`: data record download
6. `elink`: Entrez link
7. `egquery`: global query
8. `espell`: spelling suggestion
9. `ecitmatch`: batch citation searching in PubMed

See "Entrez Programming Utilities Help"
(https://www.ncbi.nlm.nih.gov/books/NBK25501/) for more details. Especially,
["E-utilities Quick Start"](https://www.ncbi.nlm.nih.gov/books/NBK25500/) is a
good starting point and ["A General Introduction to the
E-utilities"](https://www.ncbi.nlm.nih.gov/books/NBK25497/) will serve useful
information about its concepts and functions. The implemented APIs are based on
the manual on January 23, 2015.
"""
module EUtils

export
    einfo,
    esearch,
    epost,
    esummary,
    efetch,
    elink,
    egquery,
    espell,
    ecitmatch

import EzXML
import JSON
import Requests

const baseURL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"


# APIs of E-utilities
# -------------------

"""
    einfo(ctx=Dict(); params...)

Retrieve a list of databases or statistics for a database.

Parameters: db, version, retmode.
"""
function einfo(ctx::Associative=empty_context(); params...)
    params = process_parameters(params, ctx)
    return Requests.get(string(baseURL, "einfo.fcgi"), query=params)
end

"""
    esearch(ctx=Dict(); params...)

Retrieve a list of UIDs matching a text query.

Parameters: db, term, usehistory, WebEnv, query_key, retstart, retmax, rettype,
retmode, sort, field, datetype, reldate, mindate, maxdate.
"""
function esearch(ctx::Associative=empty_context(); params...)
    params = process_parameters(params, ctx)
    res = Requests.get(string(baseURL, "esearch.fcgi"), query=params)
    if get(params, :usehistory, "") == "y"
        set_context!(ctx, res)
    end
    return res
end

"""
    epost(ctx=Dict(); params...)

Upload or append a list of UIDs to the Entrez History server.

Parameters: db, id, WebEnv.
"""
function epost(ctx::Associative=empty_context(); params...)
    params = process_parameters(params, ctx)
    res = Requests.post(string(baseURL, "epost.fcgi"), data=params)
    set_context!(ctx, res)
    return res
end

"""
    esummary(ctx=Dict(); params...)

Retrieve document summaries for a list of UIDs.

Parameters: db, id, query_key, WebEnv, retstart, retmax, retmode, version.
"""
function esummary(ctx::Associative=empty_context(); params...)
    params = process_parameters(params, ctx)
    return Requests.post(string(baseURL, "esummary.fcgi"), data=params)
end

"""
    efetch(ctx=Dict(); params...)

Retrieve formatted data records for a list of UIDs.

Parameters: db, id, query_key, WebEnv, retmode, rettype, retstart, retmax,
strand, seq_start, seq_stop, complexity.
"""
function efetch(ctx::Associative=empty_context(); params...)
    params = process_parameters(params, ctx)
    return Requests.post(string(baseURL, "efetch.fcgi"), data=params)
end

"""
    elink(ctx=Dict(); params...)

Retrieve UIDs linked to an input set of UIDs.

Parameters: db, dbfrom, cmd, id, query_key, WebEnv, linkname, term, holding,
datetype, reldate, mindate, maxdate.
"""
function elink(ctx::Associative=empty_context(); params...)
    params = process_parameters(params, ctx)
    return Requests.post(string(baseURL, "elink.fcgi"), data=params)
end

"""
    egquery(ctx=Dict(); params...)

Retrieve the number of available records in all databases by a text query.

Parameters: term.
"""
function egquery(ctx::Associative=empty_context(); params...)
    params = process_parameters(params, ctx)
    return Requests.get(string(baseURL, "egquery.fcgi"), query=params)
end

"""
    espell(ctx=Dict(); params...)

Retrieve spelling suggestions.

Parameters: db, term.
"""
function espell(ctx::Associative=empty_context(); params...)
    params = process_parameters(params, ctx)
    return Requests.get(string(baseURL, "espell.fcgi"), query=params)
end

"""
    ecitmatch(ctx=Dict(); params...)

Retrieve PubMed IDs that correspond to a set of input citation strings.

Parameters: db, rettype, bdata.
"""
function ecitmatch(ctx::Associative=empty_context(); params...)
    params = process_parameters(params, ctx)
    return Requests.get(string(baseURL, "ecitmatch.cgi"), query=params)
end

# Create an empty context.
function empty_context()
    return Dict{Symbol,Any}()
end

# Set :WebEnv and :query_key values to the context `ctx` from `res`.
function set_context!(ctx, res)
    if res.status != 200
        return ctx
    end

    # extract WebEnv and query_key from the response
    contenttype = res.headers["Content-Type"]
    data = String(res.data)
    if startswith(contenttype, "text/xml")
        doc = EzXML.parsexml(data)
        ctx[:WebEnv] = EzXML.content(findfirst(doc, "//WebEnv"))
        ctx[:query_key] = EzXML.content(findfirst(doc, "//QueryKey"))
    elseif startswith(contenttype, "application/json")
        dict = JSON.parse(data)
        ctx[:WebEnv] = dict["esearchresult"]["webenv"]
        ctx[:query_key] = dict["esearchresult"]["querykey"]
    end

    return ctx
end

# Process query parameters.
function process_parameters(params, ctx)
    # merge context `ctx` into `params`
    params = merge(ctx, Dict(params))

    # flatten a set of IDs into a comma-separated string
    if haskey(params, :id)
        ids = params[:id]
        if isa(ids, AbstractString)
            ids = [ids]
        end
        params[:id] = join([string(id) for id in ids], ',')
    end

    # normalize the usehistory parameter
    if haskey(params, :usehistory) && isa(params[:usehistory], Bool)
        if params[:usehistory]
            params[:usehistory] = "y"
        else
            delete!(params, :usehistory)
        end
    end

    # stringify all values
    for (key, val) in params
        params[key] = string(val)
    end

    return params
end

end
