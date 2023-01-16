// Abstract FluxMethod Class
#ifndef _FLUX_H_
#define _FLUX_H_
#include <iostream>
#include <string>
#include "Factory.H"
#include "CAMReX.H"

class CAMReXmp;

class Flux {
public:
  Flux();
  virtual void define(std::string reconstraction) = 0;
  virtual void hyperbolicUpdate(CAMReX* AMRLevel, Array<MultiFab, AMREX_SPACEDIM>& flux_arr, Array<MultiFab, AMREX_SPACEDIM>& fluxR_arr, const CellArray& U, CellArray& U1, ParameterStruct& parameters, const Real dt, const Real *dx, const Real *prob_lo, Problem *problem, const Real time) const = 0;
  virtual ~Flux() = default;

protected:
  std::string name;
  double balance;

private:   
  static constexpr const char *def_name = "Unnamed Account";
  static constexpr double def_balance = 0.0;
};
#endif
