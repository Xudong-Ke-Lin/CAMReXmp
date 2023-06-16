model=$(cd ../../ && echo ${PWD##*/})

cfl=0.4
plot_int=2500

plot_file="/data/atalanta/xk215/PHD/Harris/Harris_RK2TVDVanLeer_c10_512x256_plt"
#mpiexec -n 12 /data/tycho/xk215/MPHIL_thesis/make_models/$model/./CAMREXmp2d.gnu.MPI.ex inputs_Harris adv.cfl=$cfl amr.plot_int=$plot_int amr.plot_file=$plot_file amr.n_cell="512 256 0" num.fluid=2 num.Maxwell=2 c=10.0 lambda_d=0.1

ls -1v $plot_file*/Header | tee $plot_file".visit"

plot_file="/data/atalanta/xk215/PHD/Harris/Harris_RK2WENOCHflat_c10_512x256_plt"
#mpiexec -n 12 /data/tycho/xk215/MPHIL_thesis/make_models/$model/./CAMREXmp2d.gnu.MPI.ex inputs_Harris adv.cfl=$cfl amr.plot_int=$plot_int amr.plot_file=$plot_file amr.n_cell="512 256 0" num.fluid=3 num.Maxwell=3 c=10.0 lambda_d=0.1

ls -1v $plot_file*/Header | tee $plot_file".visit"

########################################################################################################################################################################################################################
# c=20
########################################################################################################################################################################################################################

plot_int=5000

plot_file="/data/atalanta/xk215/PHD/Harris/Harris_RK2TVDVanLeer_c20_512x256_plt"
#mpiexec -n 12 /data/tycho/xk215/MPHIL_thesis/make_models/$model/./CAMREXmp2d.gnu.MPI.ex inputs_Harris adv.cfl=$cfl amr.plot_int=$plot_int amr.plot_file=$plot_file amr.n_cell="512 256 0" num.fluid=2 num.Maxwell=2 c=20.0 lambda_d=0.05

ls -1v $plot_file*/Header | tee $plot_file".visit"

plot_file="/data/atalanta/xk215/PHD/Harris/Harris_RK2WENOCHflat_c20_512x256_plt"
#mpiexec -n 12 /data/tycho/xk215/MPHIL_thesis/make_models/$model/./CAMREXmp2d.gnu.MPI.ex inputs_Harris adv.cfl=$cfl amr.plot_int=$plot_int amr.plot_file=$plot_file amr.n_cell="512 256 0" num.fluid=3 num.Maxwell=3 c=20.0 lambda_d=0.05

ls -1v $plot_file*/Header | tee $plot_file".visit"

########################################################################################################################################################################################################################
# Very high resolution
########################################################################################################################################################################################################################

plot_int=20000

plot_file="/data/atalanta/xk215/PHD/Harris/HarrisTest_RK2TVDVanLeer_c20_2048x1024_plt"
#mpiexec -n 16 /data/tycho/xk215/MPHIL_thesis/make_models/$model/./CAMREXmp2d.gnu.MPI.ex inputs_Harris adv.cfl=$cfl amr.plot_int=$plot_int amr.plot_file=$plot_file amr.n_cell="2048 1024 0" num.fluid=2 num.Maxwell=2 c=20.0 lambda_d=0.05 amr.max_grid_size=32

ls -1v $plot_file*/Header | tee $plot_file".visit"

plot_file="/data/atalanta/xk215/PHD/Harris/HarrisTest_RTF0.01_RK2TVDVanLeer_c20_2048x1024_plt"
mpiexec -n 16 /data/tycho/xk215/MPHIL_thesis/make_models/$model/./CAMREXmp2d.gnu.MPI.ex inputs_Harris adv.cfl=$cfl amr.plot_int=$plot_int amr.plot_file=$plot_file amr.n_cell="2048 1024 0" num.fluid=2 num.Maxwell=2 c=20.0 lambda_d=0.05 amr.max_grid_size=32

ls -1v $plot_file*/Header | tee $plot_file".visit"
