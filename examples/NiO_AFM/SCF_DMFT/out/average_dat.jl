using Printf

function read_datfile(filename)
 
    rowN = 0
    colN = 0

    firststring = ""
    open(filename, "r") do f 
        Nstring=readlines(f)
        rowN =  size(Nstring)[1]-1
        colN = size(split(Nstring[2]))[1]-1
        firststring = Nstring[1];
    end
       
    wlist = zeros(rowN)
    data_ = zeros(rowN,colN)
            
    #println(firststring);
    
    open(filename, "r") do f 
        Nstring=readlines(f)
        

        for l = 1:rowN
            wlist[l] = parse(Float64, split(Nstring[l+1])[1])
            for ll =1:colN
                data_[l,ll] = parse(Float64, split(Nstring[l+1])[ll+1])
            end
        end
        
    end
    
    return firststring, wlist, data_
end



function print_datfile(wlist, avg_data,pre_tag,suf_tag, firststring)

    filename = string(pre_tag,"avg",suf_tag)
    open(filename, "w") do io
        print(io, string("#",firststring))
        @printf(io,"\n#\n")

        for w = 1:size(wlist)[1]
            @printf(io, " %20.14f",wlist[w])
            for c = 1:size(avg_data)[2]
                @printf(io, " %20.14f", avg_data[w,c])
            end
            @printf(io,"\n")
        end
        
        
        
    end
    
    
end

println("Enter a file name to average the data (format ex. sig_smth_[1-5].dat) :")
avg_filename = readline(stdin)

rang = parse(Int64, SubString.(avg_filename, findall(r"[0-9]+",avg_filename))[1]):parse(Int64, SubString.(avg_filename, findall(r"[0-9]+",avg_filename))[2])
a=findall(r"[0-9]+",avg_filename)[1][1]
b=findall(r"[0-9]+",avg_filename)[2][end]
pre_tag = avg_filename[1:a-2]
suf_tag = avg_filename[b+2:end]


first_file = string(pre_tag,rang[1],suf_tag)


firststring, wlist, data_ = read_datfile(first_file)

avg_data = zeros(size(data_))

for i in rang
    global avg_data
    file_n = string(pre_tag,i,suf_tag)
    firststring, wlist, data_ = read_datfile(file_n);
    avg_data+=data_./size(rang)[1]
    
end

print_datfile(wlist, avg_data,pre_tag,suf_tag, firststring)

outputfilename = string(pre_tag,"avg",suf_tag)
println("\n================================================================================")
println("File '",outputfilename,"' is made successfully with files [ ",pre_tag,string(rang),suf_tag," ] !")
println("================================================================================\n")



