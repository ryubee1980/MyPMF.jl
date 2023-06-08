"""
Module for calculating potential of mean force (PMF) on the basis of Jarzynski equality. Neccessary input is a set of sample trajactory data, specifically, the pulled coodinate and the external work as functions of time. 
"""
module MyPMF
#Copyright (c) 2021 Ryuichi Okamoto <ryubee@gmail.com>
#License: https://opensource.org/licenses/MIT


"""
    pmf_HS(traj::Array{Float64,3}, ks::Float64, v::Float64, T::Float64 ; L=1000 ::Int64, energy_unit="kcal/mol" ::String)

    traj[:,:,:] is a 3-dimensional Array of the size K x T x 3, where K is the number of sample trajectories and T is the number of time slices. The time slices must be the same for all samples, traj[i,:,1]=traj[j,:,1] for all i and j.

    ks is the spring constant of biasing harmonic potential.
    v is the (linear) velocity of the biasing potential.
    T is the absolute temperature in units of Kelvin (K).
    L is the number of output data points.

    The energy_unit must be either "kcal/mol" or "kJ/mol".
    The units of length (l), time (T) can be anything, but they should consistently be used for all the variables and parameters. For example, if we set l=nm, T=ps, and E=kJ/mol, then the units of velocity and the spring constant must be [v]=nm/ps and k=kJ/mol/nm^2.

    If J_est=1, it will also calculate the Jarzynski estimation (default value is 0). That is, e^{-F/kT} is simply estimated as the arithmetic mean of e^{-w/kT}. Note that the free energy F estimated by this scheme is equal to the PMF only if the spring constant is large enough (stiff-spring limit).
"""
function pmf_HS(traj,ks,v,T; L=500,energy_unit="kcal/mol", J_est=0)
    
    if energy_unit=="kcal/mol"
        kT=T*0.593/298
    elseif energy_unit=="kJ/mol"
        kT=T*4.1868*0.593/298
    else
        throw(DomainError(energy_unit, "Energy_unit must be either kcal/mol or kJ/mol"))
    end


    # time slices. t[1] must be 0.
    t=@view traj[1,:,1]

    # The value of the coodinate coupled to the spring.
    z=@view traj[:,:,2] 

    # The (accumulated) external work during the time interval (0,t).
    w=@view traj[:,:,3] 

    # position of the potential center at time t.
    vt=@. t*v

    dl=vt[end]/L

    T=length(t)
    K=length(traj[:,1,1])

    
    window=[l*dl for l in 1:L]

    zz=@. window - dl*0.5
    
    h=zeros(Float64,T,L)
    for l in 1:L
        for i in 1:T
            for k in 1:K
                if window[l]-dl < z[k,i]-z[k,1] < window[l]
                    h[i,l]+=exp(-w[k,i]/kT)
                end
            end
        end
    end

    @. h=h/K

    
    u=Array{Float64}(undef,(T,L))

    for l in 1:L
        for i in 1:T
            
            u[i,l]=0.5*ks*(zz[l]-vt[i])^2
        end
    end
    
    eta=zeros(Float64,T)

    for i in 1:T
        for k in 1:K
            eta[i]+=exp(-w[k,i]/kT)
        end
    end
    @. eta=eta/K

    G=similar(zz)
    for l in 1:L
        temp1=0
        temp2=0
        for i in 1:T
            temp1+=h[i,l]/eta[i]
            temp2+=exp(-u[i,l]/kT)/eta[i]
        end
        G[l]=-kT*log(temp1/temp2)
    end

    Gb=G[1]

    if J_est==1
        tmp=zeros(length(vt))
        for k in 1:K
            @. tmp += exp(-(w[k,:])/kT)
           
        end
    
        pmfJ=similar(tmp)

        @. tmp = tmp/K
        @. pmfJ = -kT*log(tmp)

        return zz,G .-Gb, vt,pmfJ
    else
        return zz,G .- Gb
    end

end

"""
Essentially the same as pmf_HS, but the only difference is that pmf_HS_norm gives the PMF divided by the thermal energy k_BT.
"""
function pmf_HS_norm(traj,ks,v,T; L=500,energy_unit="kcal/mol", J_est=0)
    
    if energy_unit=="kcal/mol"
        kT=T*0.593/298
    elseif energy_unit=="kJ/mol"
        kT=T*4.1868*0.593/298
    else
        throw(DomainError(energy_unit, "Energy_unit must be either kcal/mol or kJ/mol"))
    end


    # time slices. t[1] must be 0.
    t=@view traj[1,:,1]

    # The value of the coodinate coupled to the spring.
    z=@view traj[:,:,2] 

    # The (accumulated) external work during the time interval (0,t).
    w=@view traj[:,:,3] 

    # position of the potential center at time t.
    vt=@. t*v

    dl=vt[end]/L

    T=length(t)
    K=length(traj[:,1,1])

    
    window=[l*dl for l in 1:L]

    zz=@. window - dl*0.5
    
    h=zeros(Float64,T,L)
    for l in 1:L
        for i in 1:T
            for k in 1:K
                if window[l]-dl < z[k,i]-z[k,1] < window[l]
                    h[i,l]+=exp(-w[k,i]/kT)
                end
            end
        end
    end

    @. h=h/K

    
    u=Array{Float64}(undef,(T,L))

    for l in 1:L
        for i in 1:T
            
            u[i,l]=0.5*ks*(zz[l]-vt[i])^2
        end
    end
    
    eta=zeros(Float64,T)

    for i in 1:T
        for k in 1:K
            eta[i]+=exp(-w[k,i]/kT)
        end
    end
    @. eta=eta/K

    G=similar(zz)
    for l in 1:L
        temp1=0
        temp2=0
        for i in 1:T
            temp1+=h[i,l]/eta[i]
            temp2+=exp(-u[i,l]/kT)/eta[i]
        end
        G[l]=-log(temp1/temp2)
    end

    Gb=G[1]

    if J_est==1
        tmp=zeros(length(vt))
        for k in 1:K
            @. tmp += exp(-(w[k,:])/kT)
           
        end
    
        pmfJ=similar(tmp)

        @. tmp = tmp/K
        @. pmfJ = -log(tmp)

        return zz,G .-Gb, vt,pmfJ
    else
        return zz,G .- Gb
    end

end


end
