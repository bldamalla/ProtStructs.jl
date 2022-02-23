## Parser.jl --- functions for parsing PDB files

function Base.read(filename::AbstractString, ::Type{StructureFrame}; kwargs...)
    @assert endswith(filename, ".pdb") "only `.pdb` files are supported for read"
    open(filename) do io
        read(io, StructureFrame; kwargs...)
    end
end

function Base.read(io::IO, ::Type{StructureFrame}; ensemble=false)
    ## for each in IO
    # check if it has MODEL record and create frame components
    # populate the frame components accordingly
    # read only the first model if ensemble is false
    ensemble && throw(ArgumentError("ensemble = true is not yet supported"))

    ## create frame components: step = 0, at_pos, at_list, res_list, conn
    at_pos   = SVector{3,Float64}[]
    at_list  = JAtom{Float64}[]
    res_list = JResidue[]
    conn = JConnectivity(NTuple{2,UInt}[], NTuple{3,UInt}[], 
                         NTuple{4,UInt}[], NTuple{4,UInt}[])

    # read the records
    supposed_to_be_ensemble = false
    finished_with_atoms = false
    ## residue markers
    curr_residue = :o; curr_chain = :k; curr_seqnum = 0; curr_pdb = false
    res_atom_dict = Dict{Symbol,UInt}()
    for (i, line) in enumerate(eachline(io))
        supposed_to_be_ensemble |= startswith(line, "MODEL")

        ## if starts with any of HETATM or ATOM create a JAtom
        ## at the same time take note of residues; when traversed, create a residue
        if startswith(line, "ATOM") || startswith(line, "HETATM")
            jatom, pos, res_, chain_, seqnum_, pdb_ = parseatomline(line)

            ## case for "disordered atom"; use only default one
            isnothing(jatom) && continue
            push!(at_pos, pos); push!(at_list, jatom)

            ## check both to make sure
            if res_ != curr_residue || chain_ != curr_chain || curr_seqnum != seqnum_
                # do the necessary updates
                if !isempty(res_atom_dict)
                    ch_str = curr_chain |> String
                    jresidue = JResidue(curr_residue, ch_str, res_atom_dict, curr_pdb)
                    push!(res_list, jresidue)
                end
                curr_residue = res_; curr_chain = chain_
                curr_seqnum = seqnum_; curr_pdb = pdb_
                res_atom_dict = Dict{Symbol,UInt}() ## allocate new dict
            end

            ## add atom to its residue
            res_atom_dict[jatom.name] = length(at_pos)
        end

        ## small hack: in reaching connectivity, push last recorded residue
        if !finished_with_atoms && startswith(line, "CONECT")
            ch_str = curr_chain |> String
            jresidue = JResidue(curr_residue, ch_str, res_atom_dict, curr_pdb)
            push!(res_list, jresidue)
            finished_with_atoms = true
        end

        ## TODO: HANDLE CONNECTIVITY

        if supposed_to_be_ensemble && !ensemble
            ## for now terminate early if found endmdl
            startswith(line, "ENDMDL") && break
        end
    end

    return StructureFrame(0, at_pos, at_list, res_list, conn)
end

function parseatomline(line)
    ## create an atom
    ## also return information about the chain and the residue
    (line[17] != 'A' && line[17] != ' ') && return ntuple(i->nothing, 5)

    pdb_ = startswith(line, "ATOM")
    name_ = line[13:16] |> strip |> q->replace(q, '\''=>'p') |> Symbol
    
    res_name = line[18:20] |> Symbol
    chainid = line[22] |> Symbol
    seqnum = line[23:26] |> strip |> q->parse(Int, q)

    pos_str = line[31:38], line[39:46], line[47:54]
    pos_vec = map(pos_str) do xistr
        ## try to parse
        xistr |> strip |> q->parse(Float64, q)
    end |> SVector{3,Float64}

    ## for now initialize mass and charge to 0; they're not needed anyway
    ## TODO: as for type; it's not needed but should be inferred directly
    jatom = JAtom(name_, "", zero(Float64), zero(Float64))
    jatom, pos_vec, res_name, chainid, seqnum, pdb_
end

