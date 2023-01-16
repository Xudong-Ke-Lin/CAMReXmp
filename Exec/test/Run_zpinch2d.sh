model=$(cd ../../ && echo ${PWD##*/})

stop_time=1.0
plot_file="/data/melete/xk215/MPHIL_thesis/"$model"/plt_zpinch_hyp_anex_correctdir_"
mpiexec -n 12 /data/tycho/xk215/MPHIL_thesis/make_models/$model/./CAMREXmp2d.gnu.MPI.ex inputs_zpinch2d amr.plot_file=$plot_file stop_time=$stop_time num.RK=1

#~/AMReX/amrex/Tools/Plotfile/./fextract.gnu.ex -d 1 $plot_file*/

#plot_file="/data/melete/xk215/MPHIL_thesis/"$model"/plt_zpinch_im_anex_correctdir_"
#mpiexec -n 16 /data/tycho/xk215/MPHIL_thesis/make_models/$model/./CAMREXmp2d.gnu.MPI.ex inputs_zpinch2d amr.plot_file=$plot_file stop_time=$stop_time num.RK=1

#~/AMReX/amrex/Tools/Plotfile/./fextract.gnu.ex -d 1 $plot_file


#for f in "/local/data/public/xk215/MPHIL_thesis/"$model"/plt_cyl_explosion"*; do ~/AMReX/amrex/Tools/Plotfile/./fextract.gnu.ex $f; done

