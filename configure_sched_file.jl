function read_sched_file(file)
    well_dict = Dict() #initialize dictionary
    open(file) do f
        inWELSPECS = false     #initialize variable
        for line in eachline(file)
            #println(line)
            if occursin("WELSPECS", line) #look for well's coordinates
                inWELSPECS = true
            end
            if (inWELSPECS == true) & (line == "/") #if the WELSPECS section is over
                inWELSPECS = false
            end
            if inWELSPECS == true & !occursin("WELSPECS", line)
                aux = split(line)
                push!(well_dict, aux[1] => [parse(Float64,aux[3]), parse(Float64,aux[4]), aux[5]])
            end
        end
    end
    return well_dict
end

function write_sched_file(file, well_dict)
    #Important: the number and the name of the wells is not changed

    (tmppath, tmpio) = mktemp()
    open(file) do io
        inWELSPECS = false     #initialize variable
        inCOMPDAT = false      #initialize variable
        for line in eachline(io, keep=true) # keep = true so the new line isn't chomped
            #println(line)
            #Check if line is inside the WELSPECS section----------------------
            if occursin("WELSPECS", line) #look for well's coordinates
                inWELSPECS = true
            end
            if (inWELSPECS == true) & (line == "/\n") #if the WELSPECS section is over
                inWELSPECS = false
            end

            #Check if line is inside the COMPDAT section----------------------
            if occursin("COMPDAT", line) #look for well's coordinates
                inCOMPDAT = true
            end
            if (inCOMPDAT == true) & (line == "/\n") #if the WELSPECS section is over
                inCOMPDAT = false
            end


            #Make modifications in the WELSPECS section----------------------
            if inWELSPECS == true & !occursin("WELSPECS", line)
                aux = split(line)
                well = aux[1]
                if haskey(well_dict,well) #sanity check
                    aux[3] = string(floor(Int,well_dict[well][1])) #if integer is needed: string(floor(Int,well_dict[well][1]))
                    aux[4] = string(floor(Int,well_dict[well][2]))
                    aux[5] = well_dict[well][3]
                    push!(aux, "\n")
                    line = join(aux, " ")
                end #of if well_dict[well]
            end

            #Make modifications in the COMPDAT section----------------------
            if inCOMPDAT == true & !occursin("COMPDAT", line)
                aux = split(line)
                well = aux[1]
                if haskey(well_dict,well) #sanity check
                    aux[2] = string(floor(Int,well_dict[well][1])) #if integer is needed: string(floor(Int,well_dict[well][1]))
                    aux[3] = string(floor(Int,well_dict[well][2]))
                    aux[4] = well_dict[well][3]
                    aux[5] = well_dict[well][3]
                    push!(aux, "\n")
                    line = join(aux, " ")
                end #of if well_dict[well]
            end
            write(tmpio, line)
        end
    end
    close(tmpio)
    mv(tmppath, file, force=true)
end


cd(dirname(@__FILE__)) #change location to current directory
well_dict = read_sched_file("base.sched")
well_dict["I02"][1] = 100
well_dict["P03"][2] = 200
write_sched_file("base.sched", well_dict)
