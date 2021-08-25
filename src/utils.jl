module Utils

export pwf

function pwf(s)#print with fieldnames
    for f in fieldnames(typeof(s))
        println("$f: $(getfield(s,f))")
    end
end

end