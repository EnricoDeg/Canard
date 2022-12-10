#! /bin/bash
#____________________________________________________________________________________________________
#

while getopts n:o:e: argv
do
    case "${argv}" in
        n) mpi_total_procs=${OPTARG};;
        o) io_tasks=${OPTARG};;
        e) executable=${OPTARG};;
    esac
done

set -eu

lrank=$SLURM_LOCALID

export OMPI_MCA_pml=ucx
export OMPI_MCA_btl="^vader,tcp,openib,smcuda"

# need to check in run script that the variables make sense and are
# exported!

(( compute_tasks = mpi_total_procs - io_tasks ))

if (( SLURM_PROCID < compute_tasks ))
then

    echo Compute process $SLURM_LOCALID on $(hostname)

    numanode=(2-3 0-1 6-7 4-5)
    gpus=(0 1 2 3)
    nics=(mlx5_0:1 mlx5_0:1 mlx5_1:1 mlx5_1:1)
    reorder=(0 1 2 3)

    nic_reorder=(${nics[${reorder[0]}]}
                 ${nics[${reorder[1]}]}
                 ${nics[${reorder[2]}]}
                 ${nics[${reorder[3]}]})
    numanode_reorder=(${numanode[${reorder[0]}]}
                      ${numanode[${reorder[1]}]}
                      ${numanode[${reorder[2]}]}
                      ${numanode[${reorder[3]}]})

    export UCX_NET_DEVICES=${nic_reorder[lrank]}
    export CUDA_VISIBLE_DEVICES=${gpus[${reorder[lrank]}]}

    export UCX_RNDV_SCHEME=put_zcopy
    export UCX_RNDV_THRESH=16384

    export UCX_IB_GPU_DIRECT_RDMA=yes

    export UCX_TLS=cma,rc,mm,cuda_ipc,cuda_copy,gdr_copy
    export UCX_MEMTYPE_CACHE=n

else

    echo IO process $SLURM_LOCALID on $(hostname)

    numanode=(2-3 0-1 6-7 4-5)
    nics=(mlx5_0:1 mlx5_0:1 mlx5_1:1 mlx5_1:1)    
    reorder=(0 1 2 3)
 
    nic_reorder=(${nics[${reorder[0]}]}
                 ${nics[${reorder[1]}]}
                 ${nics[${reorder[2]}]}
                 ${nics[${reorder[3]}]}) 

    numanode_reorder=(${numanode[${reorder[0]}]}
                      ${numanode[${reorder[1]}]}
                      ${numanode[${reorder[2]}]}
                      ${numanode[${reorder[3]}]}) 
    
    export UCX_NET_DEVICES=${nic_reorder[lrank]}

    export UCX_RNDV_SCHEME=put_zcopy
    export UCX_RNDV_THRESH=16384

    export UCX_IB_GPU_DIRECT_RDMA=yes

    export UCX_TLS=cma,rc,mm,cuda_ipc,cuda_copy,gdr_copy    
    export UCX_MEMTYPE_CACHE=n
    
fi

echo "numactl --cpunodebind=${numanode_reorder[$lrank]} --membind=${numanode_reorder[$lrank]} $executable"
numactl --cpunodebind=${numanode_reorder[$lrank]} --membind=${numanode_reorder[$lrank]} $executable

