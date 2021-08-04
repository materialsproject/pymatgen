# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific aliases and functions
alias sqme="squeue -u alex_epstein"
alias cdscratch="cd /global/scratch/alex_epstein"
alias loadintel="module load intel/2018.5.274.par"
alias interactive="srun -A lr_mp -p cf1 --qos condo_mp_cf1 -t 2:0:0 --pty bash -i"
alias cdconfig='cd /global/home/users/alex_epstein/.conda/envs/cms/config/'

#User specific environment variables
export FW_CONFIG_FILE='/global/home/users/alex_epstein/.conda/envs/cms/config/FW_config.yaml'
export MODULEPATH=${MODULEPATH}:/clusterfs/mp/modules

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/global/software/sl-7.x86_64/modules/langs/python/3.6/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/global/software/sl-7.x86_64/modules/langs/python/3.6/etc/profile.d/conda.sh" ]; then
        . "/global/software/sl-7.x86_64/modules/langs/python/3.6/etc/profile.d/conda.sh"
    else
        export PATH="/global/software/sl-7.x86_64/modules/langs/python/3.6/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

