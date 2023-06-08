using MyPMF
using Test

@testset "MyPMF.jl" begin
    # Write your tests here.

    K=20
    t=0.0:2.0:36000.0
    z=Array{Float64}(undef,K,length(t))
    w=similar(z)
    for k in 1:K
        z[k,:]=@. t*0.0001+0.1*(rand()-0.5)
        w[k,:]=@. (t/10000)^3-(t/6000)^2+0.5*(rand()-0.5)
    end

    dat=Array{Float64}(undef,K,length(t),3)

    @. dat[:,:,2]=z
    @. dat[:,:,3]=w

    for k in 1:K
        @. dat[k,:,1]=t
    end

    MyPMF.pmf_HS(dat,1000,0.0001,300; L=500,energy_unit="kcal/mol", J_est=0)
    MyPMF.pmf_HS(dat,1000,0.0001,300; L=500,energy_unit="kcal/mol", J_est=1)

    MyPMF.pmf_HS_norm(dat,1000,0.0001,300; L=500,energy_unit="kcal/mol", J_est=0)
    MyPMF.pmf_HS_norm(dat,1000,0.0001,300; L=500,energy_unit="kcal/mol", J_est=1)



end
