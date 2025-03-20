using JWAS,DataFrames,CSV,Statistics

global method = "BayesC"
global chain_length = 5000
global burnin = 1000

# Directly parse arguments from ARGS
for (i, arg) in enumerate(ARGS)
    if arg == "-method"
        global method = ARGS[i + 1]
    elseif arg == "-chain_length"
        global chain_length = parse(Int, ARGS[i + 1])
    elseif arg == "-burnin"
        global burnin = parse(Int, ARGS[i + 1])
    end
end

cd("/home/fotis/phd/analysis/BLUP/Sweden/single_step")
print(pwd())

pedigree = get_pedigree("SwePed.txt", separator="\t", header=true)

phenos = CSV.read("../sperm_pheno_final2024.csv", DataFrame, missingstring="NA", types=Dict(1=>String,2=>String,3=>String))

println(first(phenos))

genotypes = get_genotypes("genofile.csv", method=method, estimatePi = true, MAF=0.05)

phenos.lgDen = log.(phenos.Sperm_concentration .+ 1)

phenos[!,:VCL]= convert(Vector{Union{Missing,Float64}}, phenos[!,:VCL])
phenos[!,:lgDen]= convert(Vector{Union{Missing,Float64}}, phenos[!,:lgDen])
phenos[!,:STR]= convert(Vector{Union{Missing,Float64}}, phenos[!,:rf_EBV])

model_equation = "Sperm_concentration = intercept  + Year + Day + Id + genotypes
                    VCL = intercept  + Year + Day + Id + genotypes
                    STR = intercept + Day + Id + genotypes"

model = build_model(model_equation)

set_random(model, "Id", pedigree)
set_covariate(model, "Day")

out=runMCMC(model, phenos, chain_length=chain_length, burnin=burnin, single_step_analysis=true, 
            pedigree=pedigree, double_precision=true)
