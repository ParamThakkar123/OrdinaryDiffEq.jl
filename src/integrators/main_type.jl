type DEOptions{uEltype,uEltypeNoUnits,tTypeNoUnits,tType,F2,F3,F4,F5}
  maxiters::Int
  timeseries_steps::Int
  save_timeseries::Bool
  adaptive::Bool
  abstol::uEltype
  reltol::uEltypeNoUnits
  gamma::tTypeNoUnits
  qmax::tTypeNoUnits
  qmin::tTypeNoUnits
  dtmax::tType
  dtmin::tType
  internalnorm::F2
  progress::Bool
  progress_steps::Int
  progress_name::String
  progress_message::F5
  beta1::tTypeNoUnits
  beta2::tTypeNoUnits
  qoldinit::tTypeNoUnits
  dense::Bool
  saveat::Vector{tType}
  callback::F3
  isoutofdomain::F4
  calck::Bool
end

type ODEIntegrator{algType<:OrdinaryDiffEqAlgorithm,uType<:Union{AbstractArray,Number},tType,tstopsType,tTypeNoUnits,ksEltype,SolType,rateType,F,ProgressType,CacheType,ECType,O} <: AbstractODEIntegrator
  sol::SolType
  u::uType
  k::ksEltype
  t::tType
  dt::tType
  f::F
  uprev::uType
  kprev::ksEltype
  tprev::tType
  tstops::tstopsType
  saveat::tstopsType
  adaptiveorder::Int
  order::Int
  alg::algType
  rate_prototype::rateType
  notsaveat_idxs::Vector{Int}
  calcprevs::Bool
  dtcache::tType
  dtpropose::tType
  dt_mod::tTypeNoUnits
  tdir::Int
  qminc::tTypeNoUnits
  qmaxc::tTypeNoUnits
  EEst::tTypeNoUnits
  qold::tTypeNoUnits
  iter::Int
  saveiter::Int
  saveiter_dense::Int
  prog::ProgressType
  cache::CacheType
  event_cache::ECType
  kshortsize::Int
  reeval_fsal::Bool
  advance_to_tstop::Bool
  opts::O
  fsalfirst::rateType
  fsallast::rateType

  ODEIntegrator(sol,u,k,t,dt,f,uprev,kprev,tprev,tstops,saveat,adaptiveorder,
    order,alg,rate_prototype,notsaveat_idxs,calcprevs,dtcache,dtpropose,dt_mod,tdir,qminc,
    qmaxc,EEst,qold,iter,saveiter,saveiter_dense,prog,cache,event_cache,
    kshortsize,reeval_fsal,advance_to_tstop,opts) = new(
    sol,u,k,t,dt,f,uprev,kprev,tprev,tstops,saveat,adaptiveorder,
      order,alg,rate_prototype,notsaveat_idxs,calcprevs,dtcache,dtpropose,dt_mod,tdir,qminc,
      qmaxc,EEst,qold,iter,saveiter,saveiter_dense,prog,cache,event_cache,
      kshortsize,reeval_fsal,advance_to_tstop,opts) # Leave off fsalfirst and last
end

function start(integrator::ODEIntegrator)
  initialize!(integrator,integrator.cache)
  integrator.iter
end

function next(integrator::ODEIntegrator,state)
  state += 1
  if integrator.advance_to_tstop
    while integrator.tdir*integrator.t < integrator.tdir*top(integrator.tstops)
      ode_loopheader!(integrator)
      perform_step!(integrator,integrator.cache)
      ode_loopfooter!(integrator)
    end
  else
    ode_loopheader!(integrator)
    perform_step!(integrator,integrator.cache)
    ode_loopfooter!(integrator)
  end
  integrator.t == top(integrator.tstops) && pop!(integrator.tstops)
  integrator,state
end

function done(integrator::ODEIntegrator,state)
  if integrator.iter > integrator.opts.maxiters
    warn("Interrupted. Larger maxiters is needed.")
    ode_postamble!(integrator)
    return true
  end
  if integrator.dt == zero(integrator.t)
    warn("dt == 0. Aborting")
    ode_postamble!(integrator)
    return true
  end
  if any(isnan,integrator.uprev)
    warn("NaNs detected. Aborting")
    ode_postamble!(integrator)
    return true
  end
  if isempty(integrator.tstops)
    ode_postamble!(integrator)
    return true
  end
  false
end

eltype(integrator::ODEIntegrator) = typeof(integrator)
