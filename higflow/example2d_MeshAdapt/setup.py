#!/usr/bin/env python3
import math
import re
from sys import argv
import ruamel.yaml
yaml = ruamel.yaml.YAML()

filename = 'input/load.par.contr.yaml'

# Parse arguments
args = argv[1:]
if(len(args)>0):
    if('help' in args[0]):
        print(f'Usage: \nEither (for singlephase) \n./setup [type] [model](of visc or eo - optional) [params](visc, eo or general) ...')
        print(f'or (for multiphase) \n./setup mult [type0-type1] [model0-model1](of visc or eo - optional) [params](visc, eo or general) ...')
        print('type: newt, visc, eo, eovisc (or visceo)')
        exit()

def pos_in_args(str):
    for pos, element in enumerate(args):
        if str in element:
            return pos
    else:
        return -1
 
with open(filename) as yaml_file:
    yaml_data = yaml.load(yaml_file)

################## Change the controllers and parameters based on input ##################

def set_param(yaml_dict, keys, input_str, cast2int=False):
    par_pos = pos_in_args(input_str+'=')
    if(par_pos != -1):
        current = yaml_dict
        for key in keys[:-1]:
            current = current[key]
        # Assign the value to the last key
        if cast2int:
            current[keys[-1]] = int(args[par_pos].split('=')[1])
        else:
            current[keys[-1]] = float(args[par_pos].split('=')[1])

def check_and_assign(yaml_dict, keys, value):
    current = yaml_dict
    for key in keys[:-1]:
        try:
            current = current[key]
        except:
            return False
    try:
        current[keys[-1]] = value
    except:
        return False
    return True

def check_and_compare(yaml_dict, keys, value, default):
    current = yaml_dict
    for key in keys[:-1]:
        try:
            current = current[key]
        except:
            return default
    try:
        return (current[keys[-1]] == value)
    except:
        return default
    
def get_param_or_default(yaml_dict, keys, default):
    current = yaml_dict
    for key in keys[:-1]:
        try:
            current = current[key]
        except:
            return default
    try:
        return current[keys[-1]]
    except:
        return default

if len(args) > 0:

    if args[0].strip() == 'mult':
        def found_model(model, phase):
            if model+str(phase) in args:
                return True
            model_pos = pos_in_args(model)
            if model_pos != -1:
                model_arg = args[model_pos]
                if len(model_arg.split('-')) == 2:
                    if(phase == 0):
                        return model == model_arg.split('-')[0]
                    elif(phase == 1):
                        return model == model_arg.split('-')[1]
            return False

        alphaeo_can_be_called_alpha = True
        yaml_data['flowphase'] = 'multiphase'
        if(len(args) < 2):
            flowtype0 = 'newt'
            flowtype1 = 'newt'
        else:
            split_arg1 = args[1].split('-')

            if(len(split_arg1)>1):
                flowtype0 = split_arg1[0]
                flowtype1 = split_arg1[1]
            else:
                flowtype0 = split_arg1[0]
                flowtype1 = 'newt'

        if(check_and_compare(yaml_data, ['multiphase','contr','surf_tension'], True, True)):
            set_param(yaml_data, ['multiphase', 'adimensional', 'Ca'], 'Ca')
        par_pos = pos_in_args('stension')
        if(par_pos != -1):
            yaml_data['multiphase']['contr']['surf_tension'] = True

        par_pos = pos_in_args('nostension')
        if(par_pos != -1):
            yaml_data['multiphase']['contr']['surf_tension'] = False

        set_param(yaml_data, ['multiphase', 'phase0', 'adimensional', 'rho'], 'rho0')
        set_param(yaml_data, ['multiphase', 'phase1', 'adimensional', 'rho'], 'rho1')
        set_param(yaml_data, ['multiphase', 'phase0', 'adimensional', 'mu'], 'mu0')
        set_param(yaml_data, ['multiphase', 'phase1', 'adimensional', 'mu'], 'mu1')
        

        both_visc = 'visc' in flowtype0 and 'visc' in flowtype1

            
        if 'newt' in flowtype0:
            yaml_data['multiphase']['contr']['flowtype0'] = 'newtonian'
        elif 'visc' in flowtype0:
            yaml_data['multiphase']['contr']['flowtype0'] = 'viscoelastic'

            visc_dict = yaml_data['multiphase']['viscoelastic']
            visc_dict0 = yaml_data['multiphase']['viscoelastic']['phase0']
            
            set_param(visc_dict0, ['adimensional', 'De'], 'De0')
            set_param(visc_dict0, ['adimensional', 'beta'], 'beta0')
            if(not both_visc):
                set_param(visc_dict0, ['adimensional', 'De'], 'De')
                set_param(visc_dict0, ['adimensional', 'beta'], 'beta')
            
            if (found_model('oldroyd',0) or found_model('oldroyd_b',0) or 
                (not both_visc and ('oldroyd' in args or 'oldroyd_b' in args))):
                visc_dict['contr']['model0'] = 'oldroyd_b'
            elif found_model('giesekus',0) or (not both_visc and 'giesekus' in args):
                visc_dict['contr']['model0'] = 'giesekus'
                set_param(visc_dict0, ['model_giesekus', 'alpha'], 'alpha0')
                if(not (both_visc and found_model('giesekus',1))):
                    set_param(visc_dict0, ['model_giesekus', 'alpha'], 'alpha')
            elif found_model('lptt',0) or (not both_visc and 'lptt' in args):
                visc_dict['contr']['model0'] = 'lptt'
                set_param(visc_dict0, ['model_lptt', 'epsilon'], 'eps0')
                set_param(visc_dict0, ['model_lptt', 'xi'], 'xi0')
                if(not (both_visc and (found_model('lptt',1) or found_model('gptt',1)))):
                    set_param(visc_dict0, ['model_lptt', 'epsilon'], 'eps')
                    set_param(visc_dict0, ['model_lptt', 'xi'], 'xi')
                alphaeo_can_be_called_alpha = False
            elif found_model('gptt',0) or (not both_visc and 'gptt' in args):
                visc_dict['contr']['model0'] = 'gptt'
                set_param(visc_dict0, ['model_gptt', 'epsilon'], 'eps0')
                set_param(visc_dict0, ['model_gptt', 'xi'], 'xi0')
                set_param(visc_dict0, ['model_gptt', 'alpha_gptt'], 'alpha_gptt0')
                set_param(visc_dict0, ['model_gptt', 'beta_gptt'], 'beta_gptt0')
                if(not (both_visc and (found_model('gptt',1) or found_model('lptt',1)))):
                    set_param(visc_dict0, ['model_gptt', 'epsilon'], 'eps')
                    set_param(visc_dict0, ['model_gptt', 'xi'], 'xi')
                    set_param(visc_dict0, ['model_gptt', 'alpha_gptt'], 'alpha_gptt')
                    set_param(visc_dict0, ['model_gptt', 'beta_gptt'], 'beta_gptt')
                alphaeo_can_be_called_alpha = False
            elif found_model('fene_p',0) or (not both_visc and 'fene_p' in args):
                visc_dict['contr']['model0'] = 'fene_p'
                set_param(visc_dict0, ['model_fene_p', 'L2'], 'L20')
                set_param(visc_dict0, ['model_fene_p', 'L2'], 'Lsq0')
                if(not (both_visc and (found_model('fene_p',1) or found_model('e_fene',1)))):
                    set_param(visc_dict0, ['model_fene_p', 'L2'], 'L2')
                    set_param(visc_dict0, ['model_fene_p', 'L2'], 'Lsq')
            elif found_model('e_fene',0) or (not both_visc and 'e_fene' in args):
                visc_dict['contr']['model0'] = 'e_fene'
                set_param(visc_dict0, ['model_e_fene', 'L2'], 'L20')
                set_param(visc_dict0, ['model_e_fene', 'L2'], 'Lsq0')
                set_param(visc_dict0, ['model_e_fene', 'lambda'], 'lambda0')
                set_param(visc_dict0, ['model_e_fene', 'E'], 'E0')
                if(not (both_visc and (found_model('e_fene',1) or found_model('fene_p',1)))):
                    set_param(visc_dict0, ['model_e_fene', 'L2'], 'L2')
                    set_param(visc_dict0, ['model_e_fene', 'L2'], 'Lsq')
                    set_param(visc_dict0, ['model_e_fene', 'lambda'], 'lambda')
                    set_param(visc_dict0, ['model_e_fene', 'E'], 'E')

        if 'newt' in flowtype1:
            yaml_data['multiphase']['contr']['flowtype1'] = 'newtonian'
        elif 'visc' in flowtype1:
            yaml_data['multiphase']['contr']['flowtype1'] = 'viscoelastic'
            visc_dict = yaml_data['multiphase']['viscoelastic']
            visc_dict1 = yaml_data['multiphase']['viscoelastic']['phase1']

            set_param(visc_dict1, ['adimensional', 'De'], 'De1')
            set_param(visc_dict1, ['adimensional', 'beta'], 'beta1')
            if(not both_visc):
                set_param(visc_dict1, ['adimensional', 'De'], 'De')
                set_param(visc_dict1, ['adimensional', 'beta'], 'beta')
            
            if (found_model('oldroyd',1) or found_model('oldroyd_b',1) or
                (not both_visc and ('oldroyd' in args or 'oldroyd_b' in args))):
                visc_dict['contr']['model1'] = 'oldroyd_b'
            elif found_model('giesekus',1) or (not both_visc and 'giesekus' in args):
                visc_dict['contr']['model1'] = 'giesekus'
                set_param(visc_dict1, ['model_giesekus', 'alpha'], 'alpha0')
                if(not (both_visc and found_model('giesekus',0))):
                    set_param(visc_dict1, ['model_giesekus', 'alpha'], 'alpha')
            elif found_model('lptt',1) or (not both_visc and 'lptt' in args):
                visc_dict['contr']['model1'] = 'lptt'
                set_param(visc_dict1, ['model_lptt', 'epsilon'], 'eps1')
                set_param(visc_dict1, ['model_lptt', 'xi'], 'xi1')
                if(not (both_visc and (found_model('lptt',0) or found_model('gptt',0)))):
                    set_param(visc_dict1, ['model_lptt', 'epsilon'], 'eps')
                    set_param(visc_dict1, ['model_lptt', 'xi'], 'xi')
                alphaeo_can_be_called_alpha = False
            elif found_model('gptt',1) or (not both_visc and 'gptt' in args):
                visc_dict['contr']['model1'] = 'gptt'
                set_param(visc_dict1, ['model_gptt', 'epsilon'], 'eps1')
                set_param(visc_dict1, ['model_gptt', 'xi'], 'xi1')
                set_param(visc_dict1, ['model_gptt', 'alpha_gptt'], 'alpha_gptt1')
                set_param(visc_dict1, ['model_gptt', 'beta_gptt'], 'beta_gptt1')
                if(not (both_visc and (found_model('gptt',0) or found_model('lptt',0)))):
                    set_param(visc_dict1, ['model_gptt', 'epsilon'], 'eps')
                    set_param(visc_dict1, ['model_gptt', 'xi'], 'xi')
                    set_param(visc_dict1, ['model_gptt', 'alpha_gptt'], 'alpha_gptt')
                    set_param(visc_dict1, ['model_gptt', 'beta_gptt'], 'beta_gptt')
                alphaeo_can_be_called_alpha = False
            elif found_model('fene_p',1) or (not both_visc and 'fene_p' in args):
                visc_dict['contr']['model1'] = 'fene_p'
                set_param(visc_dict1, ['model_fene_p', 'L2'], 'L21')
                set_param(visc_dict1, ['model_fene_p', 'L2'], 'Lsq1')
                if(not (both_visc and (found_model('fene_p',0) or found_model('e_fene',0)))):
                    set_param(visc_dict1, ['model_fene_p', 'L2'], 'L2')
                    set_param(visc_dict1, ['model_fene_p', 'L2'], 'Lsq')
            elif found_model('e_fene',1) or (not both_visc and 'e_fene' in args):
                visc_dict['contr']['model1'] = 'e_fene'
                set_param(visc_dict1, ['model_e_fene', 'L2'], 'L21')
                set_param(visc_dict1, ['model_e_fene', 'L2'], 'Lsq1')
                set_param(visc_dict1, ['model_e_fene', 'lambda'], 'lambda1')
                set_param(visc_dict1, ['model_e_fene', 'E'], 'E1')
                if(not (both_visc and (found_model('e_fene',0) or found_model('fene_p',0)))):
                    set_param(visc_dict1, ['model_e_fene', 'L2'], 'L2')
                    set_param(visc_dict1, ['model_e_fene', 'L2'], 'Lsq')
                    set_param(visc_dict1, ['model_e_fene', 'lambda'], 'lambda')
                    set_param(visc_dict1, ['model_e_fene', 'E'], 'E')


        both_eo = 'eo' in flowtype0 and 'eo' in flowtype1

        if 'eo' in flowtype0:
            yaml_data['multiphase']['contr']['eoflow0'] = True
            if 'visc' not in flowtype0:
                yaml_data['multiphase']['contr']['flowtype0'] = 'newtonian'
            eo_dict = yaml_data['multiphase']['electroosmotic']
            eo_dict0 = yaml_data['multiphase']['electroosmotic']['phase0']

            set_param(eo_dict0, ['adimensional', 'alpha'], 'alphaeo0')
            if(alphaeo_can_be_called_alpha):
                set_param(eo_dict0, ['adimensional', 'alpha'], 'alpha0')
            kappa_pos = pos_in_args('kappa0=')
            if(kappa_pos != -1):
                kappa = float(args[kappa_pos].split('=')[1])
                alpha = eo_dict0['adimensional']['alpha']
                eo_dict0['adimensional']['delta'] = 0.5*kappa**2/alpha  
            set_param(eo_dict0, ['adimensional', 'Ex'], 'Ex0')
            eo_dict0['adimensional']['perm'] = 1.0 # default value
            set_param(eo_dict0, ['adimensional', 'perm'], 'perm0')
            set_param(eo_dict0, ['model_pnp', 'Pe'], 'Pe0')
            if(not both_eo):
                set_param(eo_dict0, ['adimensional', 'alpha'], 'alphaeo')
                if(alphaeo_can_be_called_alpha):
                    set_param(eo_dict0, ['adimensional', 'alpha'], 'alpha')
                kappa_pos = pos_in_args('kappa=')
                if(kappa_pos != -1):
                    kappa = float(args[kappa_pos].split('=')[1])
                    alpha = eo_dict0['adimensional']['alpha']
                    eo_dict0['adimensional']['delta'] = 0.5*kappa**2/alpha 
                set_param(eo_dict0, ['adimensional', 'Ex'], 'Ex')
                set_param(eo_dict0, ['model_pnp', 'Pe'], 'Pe')
            
            if 'pnp' in args or found_model('pnp',0):
                eo_dict['contr']['model'] = 'pnp'
            if 'pb' in args or found_model('pb',0):
                eo_dict['contr']['model'] = 'pb'
            elif 'pbdh' in args or found_model('pbdh',0):
                eo_dict['contr']['model'] = 'pbdh'
            set_param(eo_dict, ['model_pnp', 'max_inner_iter'], 'inner', cast2int=True)
        else:
            check_and_assign(yaml_data, ['multiphase','contr','eoflow0'], False)


        if 'eo' in flowtype1:
            yaml_data['multiphase']['contr']['eoflow1'] = True
            if 'visc' not in flowtype1:
                yaml_data['multiphase']['contr']['flowtype1'] = 'newtonian'
            eo_dict = yaml_data['multiphase']['electroosmotic']
            eo_dict1 = yaml_data['multiphase']['electroosmotic']['phase1']

            set_param(eo_dict1, ['adimensional', 'alpha'], 'alphaeo1')
            if(alphaeo_can_be_called_alpha):
                set_param(eo_dict1, ['adimensional', 'alpha'], 'alpha1')
            kappa_pos = pos_in_args('kappa1=')
            if(kappa_pos != -1):
                kappa = float(args[kappa_pos].split('=')[1])
                alpha = eo_dict1['adimensional']['alpha']
                eo_dict1['adimensional']['delta'] = 0.5*kappa**2/alpha
            set_param(eo_dict1, ['adimensional', 'Ex'], 'Ex1')
            eo_dict1['adimensional']['perm'] = 1.0 # default value
            set_param(eo_dict1, ['adimensional', 'perm'], 'perm1')
            set_param(eo_dict1, ['model_pnp', 'Pe'], 'Pe1')
            if(not both_eo):
                set_param(eo_dict1, ['adimensional', 'alpha'], 'alphaeo')
                if(alphaeo_can_be_called_alpha):
                    set_param(eo_dict1, ['adimensional', 'alpha'], 'alpha')
                kappa_pos = pos_in_args('kappa=')
                if(kappa_pos != -1):
                    kappa = float(args[kappa_pos].split('=')[1])
                    alpha = eo_dict1['adimensional']['alpha']
                    eo_dict1['adimensional']['delta'] = 0.5*kappa**2/alpha
                set_param(eo_dict1, ['adimensional', 'Ex'], 'Ex')
                set_param(eo_dict1, ['model_pnp', 'Pe'], 'Pe')

            if 'pnp' in args or found_model('pnp',1):
                eo_dict['contr']['model'] = 'pnp'
            if 'pb' in args or found_model('pb',1):
                eo_dict['contr']['model'] = 'pb'
            elif 'pbdh' in args or found_model('pbdh',1):
                eo_dict['contr']['model'] = 'pbdh'
            set_param(eo_dict, ['model_pnp', 'max_inner_iter'], 'inner', cast2int=True)
        else:
            check_and_assign(yaml_data, ['multiphase','contr','eoflow1'], False)

    else: # singlephase
        alphaeo_can_be_called_alpha = True
        yaml_data['flowphase'] = 'singlephase'
        flowtype = args[0]
        if 'newt' in flowtype:
            yaml_data['singlephase']['contr']['flowtype'] = 'newtonian'
        elif 'visc' in flowtype:
            yaml_data['singlephase']['contr']['flowtype'] = 'viscoelastic'

            visc_dict = yaml_data['singlephase']['viscoelastic']
            set_param(visc_dict, ['adimensional', 'De'], 'De')
            set_param(visc_dict, ['adimensional', 'beta'], 'beta')

            if 'oldroyd' in args or 'oldroyd_b' in args:
                visc_dict['contr']['model'] = 'oldroyd_b'
            elif 'giesekus' in args:
                visc_dict['contr']['model'] = 'giesekus'
                set_param(visc_dict, ['model_giesekus', 'alpha'], 'alpha')
                alphaeo_can_be_called_alpha = False
            elif 'lptt' in args:
                visc_dict['contr']['model'] = 'lptt'
                set_param(visc_dict, ['model_lptt', 'epsilon'], 'eps')
                set_param(visc_dict, ['model_lptt', 'xi'], 'xi')
            elif 'gptt' in args:
                visc_dict['contr']['model'] = 'gptt'
                set_param(visc_dict, ['model_gptt', 'epsilon'], 'eps')
                set_param(visc_dict, ['model_gptt', 'xi'], 'xi')
                set_param(visc_dict, ['model_gptt', 'alpha_gptt'], 'alpha_gptt')
                set_param(visc_dict, ['model_gptt', 'beta_gptt'], 'beta_gptt')
                alphaeo_can_be_called_alpha = False
            elif 'fene_p' in args:
                visc_dict['contr']['model'] = 'fene_p'
                set_param(visc_dict, ['model_fene_p', 'L2'], 'L2')
                set_param(visc_dict, ['model_fene_p', 'L2'], 'Lsq')
            elif 'e_fene' in args:
                visc_dict['contr']['model'] = 'e_fene'
                set_param(visc_dict, ['model_e_fene', 'L2'], 'L2')
                set_param(visc_dict, ['model_fene_p', 'L2'], 'Lsq')
                set_param(visc_dict, ['model_e_fene', 'lambda'], 'lambda')
                set_param(visc_dict, ['model_e_fene', 'E'], 'E')

        if 'eo' in flowtype:
            yaml_data['singlephase']['contr']['eoflow'] = True
            if 'visc' not in flowtype:
                yaml_data['singlephase']['contr']['flowtype'] = 'newtonian'

            eo_dict = yaml_data['singlephase']['electroosmotic']
            set_param(eo_dict, ['adimensional', 'alpha'], 'alphaeo')
            if(alphaeo_can_be_called_alpha):
                set_param(eo_dict, ['adimensional', 'alpha'], 'alpha')
            kappa_pos = pos_in_args('kappa=')
            if(kappa_pos != -1):
                kappa = float(args[kappa_pos].split('=')[1])
                alpha = eo_dict['adimensional']['alpha']
                eo_dict['adimensional']['delta'] = 0.5*kappa**2/alpha
            set_param(eo_dict, ['adimensional', 'Ex'], 'Ex')
            set_param(eo_dict, ['model_pnp', 'Pe'], 'Pe')

            if 'pnp' in args:
                eo_dict['contr']['model'] = 'pnp'
            if 'pb' in args:
                eo_dict['contr']['model'] = 'pb'
            elif 'pbdh' in args:
                eo_dict['contr']['model'] = 'pbdh'
            set_param(eo_dict, ['model_pnp', 'max_inner_iter'], 'inner', cast2int=True)
        else:
            check_and_assign(yaml_data, ['singlephase','contr','eoflow'], False)

    set_param(yaml_data, ['adimensional', 'Re'], 'Re')
    if(check_and_compare(yaml_data, ['simulation_contr','gravity'], True, False)):
        set_param(yaml_data, ['adimensional', 'Fr'], 'Fr')
    set_param(yaml_data, ['simulation_par', 'numsteps'], 'numsteps', cast2int=True)
    set_param(yaml_data, ['simulation_par', 'dt'], 'dt')
    set_param(yaml_data, ['simulation_par', 'dts'], 'dts')
    set_param(yaml_data, ['simulation_par', 'dtp'], 'dtp')
    numsteps_pos = pos_in_args('numsteps=')
    if(numsteps_pos != -1):
        numsteps_ = int(args[numsteps_pos].split('=')[1])
    else:
        numsteps_ = 1
    tf_pos = pos_in_args('tf=')
    if(tf_pos != -1):
        tf = float(args[tf_pos].split('=')[1])
        dt = yaml_data['simulation_par']['dt']
        yaml_data['simulation_par']['numsteps'] = max(int(tf/dt) + 1, numsteps_)

with open(filename, 'w') as yaml_file:
    yaml.dump(yaml_data, yaml_file)


eo_any = False
if(yaml_data['flowphase'] == 'singlephase'):
    if(check_and_compare(yaml_data, ['singlephase','contr','eoflow'], True, False)):
        eo_any = True
elif(yaml_data['flowphase'] == 'multiphase'):
    if(check_and_compare(yaml_data, ['multiphase','contr','eoflow0'], True, False) or 
        check_and_compare(yaml_data, ['multiphase','contr','eoflow1'], True, False)):
        eo_any = True

######################### Change boundary conditions based on input #########################

if 'cavity' in args or 'channel' in args:
    bc_filename = 'input/load.bc.yaml'

    with open(bc_filename) as yaml_file:
        bc_yaml_data = yaml.load(yaml_file)

    if 'cavity' in args:
        bc_yaml_data['bc']['bc0']['pressure']['type'] = 'neumann'
        bc_yaml_data['bc']['bc0']['velocity_0']['type'] = 'dirichlet'
        bc_yaml_data['bc']['bc0']['velocity_1']['type'] = 'dirichlet'
        bc_yaml_data['bc']['bc1']['pressure']['type'] = 'neumann'
        bc_yaml_data['bc']['bc1']['velocity_0']['type'] = 'dirichlet'
        bc_yaml_data['bc']['bc1']['velocity_1']['type'] = 'dirichlet'
        bc_yaml_data['bc']['bc2']['pressure']['type'] = 'neumann'
        bc_yaml_data['bc']['bc2']['velocity_0']['type'] = 'dirichlet'
        bc_yaml_data['bc']['bc2']['velocity_1']['type'] = 'dirichlet'
        bc_yaml_data['bc']['bc3']['pressure']['type'] = 'neumann'
        bc_yaml_data['bc']['bc3']['velocity_0']['type'] = 'dirichlet'
        bc_yaml_data['bc']['bc3']['velocity_1']['type'] = 'dirichlet'

    
    if 'channel' in args:
        if(eo_any == False):
            bc_yaml_data['bc']['bc0']['pressure']['type'] = 'neumann'
            bc_yaml_data['bc']['bc0']['velocity_0']['type'] = 'dirichlet'
            bc_yaml_data['bc']['bc0']['velocity_1']['type'] = 'dirichlet'
        else:
            bc_yaml_data['bc']['bc0']['pressure']['type'] = 'dirichlet'
            bc_yaml_data['bc']['bc0']['velocity_0']['type'] = 'neumann'
            bc_yaml_data['bc']['bc0']['velocity_1']['type'] = 'neumann'
        bc_yaml_data['bc']['bc1']['pressure']['type'] = 'neumann'
        bc_yaml_data['bc']['bc1']['velocity_0']['type'] = 'dirichlet'
        bc_yaml_data['bc']['bc1']['velocity_1']['type'] = 'dirichlet'
        bc_yaml_data['bc']['bc2']['pressure']['type'] = 'dirichlet'
        bc_yaml_data['bc']['bc2']['velocity_0']['type'] = 'neumann'
        bc_yaml_data['bc']['bc2']['velocity_1']['type'] = 'neumann'
        bc_yaml_data['bc']['bc3']['pressure']['type'] = 'neumann'
        bc_yaml_data['bc']['bc3']['velocity_0']['type'] = 'dirichlet'
        bc_yaml_data['bc']['bc3']['velocity_1']['type'] = 'dirichlet'

    with open(bc_filename, 'w') as yaml_file:
        yaml.dump(bc_yaml_data, yaml_file)

name_eo_symmetry = ''
if(eo_any):
    defines_filename = 'ns-example-2d.h'
    psi_down_vals = {'sym': -1.0, 'asym': 1.0, 'half': 0.0}
    psi_down_value = -1.0
    for k in psi_down_vals.keys():
        if(k in args):
            if(k!= 'sym'):
                name_eo_symmetry = k+'_'
            psi_down_value = psi_down_vals[k]
            break        
    with open(defines_filename, 'r') as file:
        defines_file = file.read()
    defines_file = re.sub(
        rf'(psi_down\s+=\s+)[^\s;]+(?=;)',
        f'\\g<1>{psi_down_value}',
        defines_file,
        flags=re.DOTALL
    )
    with open(defines_filename, 'w') as file:
        file.write(defines_file)

######################### Create name of output (to return) using controllers and parameters #########################

def par2name(par_path):
    decs = 3
    if(par_path < 1e-2):
        decs = 4
    if(par_path < 1e-3):
        decs = 5
    if(par_path < 1e-4):
        decs = 6
    if(par_path < 1e-5):
        decs = 7
    return str(round(par_path, decs))


name_params = '-Re' + par2name(yaml_data['adimensional']['Re']) 
if(check_and_compare(yaml_data, ['simulation_contr','gravity'], True, False)):
    name_params += '_Fr' + par2name(yaml_data['adimensional']['Fr'])

if(yaml_data['flowphase'] == 'singlephase'):
    name_eoflow = ''
    name_eoflow_model = ''
    name_eoflow_params = ''
    name_flowtype = ''
    name_flowtype_model = ''
    name_flowtype_params = ''
    if(check_and_compare(yaml_data, ['singlephase','contr','eoflow'], True, False)):
        yaml_eo = yaml_data['singlephase']['electroosmotic']
        name_eoflow = name_eo_symmetry+'EO'
        if yaml_eo['contr']['model'] == 'pb':
            name_eoflow_model = 'pb'
        elif yaml_eo['contr']['model'] == 'pbdh':
            name_eoflow_model = 'pbdh'  
        else:
            name_eoflow_model = '' # pnp needs no description
        param_path = yaml_eo['adimensional']
        name_eoflow_params = '\u03BA' + par2name(math.sqrt(2.0*param_path['alpha']*param_path['delta']))
        name_eoflow_params += '_\u03B1' + par2name(param_path['alpha'])
        if yaml_eo['contr']['model'] == 'pnp':
            name_eoflow_params += '_Pe' + par2name(yaml_eo['model_pnp']['Pe'])
        name_eoflow_params += '_Ex' + par2name(param_path['Ex'])
    elif(yaml_data['singlephase']['contr']['flowtype'] == 'newtonian'):
        name_flowtype = 'Newt'

    if(yaml_data['singlephase']['contr']['flowtype'] == 'viscoelastic'):
        yaml_visc = yaml_data['singlephase']['viscoelastic']
        name_flowtype = 'Visc'
        name_flowtype_model = yaml_visc['contr']['model']
        name_flowtype_params = 'De' + par2name(yaml_visc['adimensional']['De']) + '_\u03B2' + par2name(yaml_visc['adimensional']['beta'])
        if(name_flowtype_model == 'giesekus'):
            name_flowtype_params += '_\u03B1' + par2name(yaml_visc['model_giesekus']['alpha'])
        elif(name_flowtype_model == 'lptt'):
            name_flowtype_params += '_\u03B5' + par2name(yaml_visc['model_lptt']['epsilon']) 
            if(yaml_visc['model_lptt']['xi'] != 0.0):
                name_flowtype_params += '_\u03BE' + par2name(yaml_visc['model_lptt']['xi'])
        elif(name_flowtype_model == 'gptt'):
            name_flowtype_params += '_\u03B5' + par2name(yaml_visc['model_gptt']['epsilon']) 
            if(yaml_visc['model_gptt']['xi'] != 0.0):
                name_flowtype_params += '_\u03BE' + par2name(yaml_visc['model_gptt']['xi'])
            if(yaml_visc['model_gptt']['alpha_gptt'] != 1.0):
                name_flowtype_params += '_\u03B1_ml' + par2name(yaml_visc['model_gptt']['alpha_gptt'])
            if(yaml_visc['model_gptt']['beta_gptt'] != 1.0):
                name_flowtype_params += '_\u03B2_ml' + par2name(yaml_visc['model_gptt']['beta_gptt'])
        elif(name_flowtype_model == 'fene_p'):
            name_flowtype_params += '_L\u00B2' + par2name(yaml_visc['model_fene_p']['L2'])
        elif(name_flowtype_model == 'e_fene'):
            name_flowtype_params += '_L\u00B2' + par2name(yaml_visc['model_e_fene']['L2'])
            name_flowtype_params += '_\u03BB' + par2name(yaml_visc['model_e_fene']['lambda']) + '_E' + par2name(yaml_visc['model_e_fene']['E'])

    if name_flowtype != '':
        name_flowtype = name_flowtype + '-'
    if name_flowtype_model != '':
        name_flowtype_model = name_flowtype_model + '-'
    if name_eoflow != '':
        name_eoflow = name_eoflow + '-'
    if name_eoflow_model != '':
        name_eoflow_model = name_eoflow_model + '-'
    if name_eoflow_params != '':
        name_eoflow_params = '-' + name_eoflow_params
    if name_flowtype_params != '':
        name_flowtype_params = '-' + name_flowtype_params

    name = f"{name_eoflow}{name_eoflow_model}{name_flowtype}{name_flowtype_model}"
    name += f"{name_eoflow_params}{name_flowtype_params}{name_params}"
    

elif(yaml_data['flowphase'] == 'multiphase'):
    name_phase = 'mult_'
    name_mult_params = ''
    name_mult0_params = ''
    name_mult1_params = ''
    if(check_and_compare(yaml_data, ['multiphase','contr','surf_tension'], True, True)):
        name_mult_params += '-Ca' + par2name(yaml_data['multiphase']['adimensional']['Ca'])

    if(get_param_or_default(yaml_data, ['multiphase', 'phase0', 'adimensional', 'rho'], 1.0) != 1.0):
        name_mult0_params += '_\u03C1' + par2name(yaml_data['multiphase']['phase0']['adimensional']['rho'])
    if(get_param_or_default(yaml_data, ['multiphase', 'phase0', 'adimensional', 'mu'], 1.0) != 1.0):
        name_mult0_params += '_\u03BC' + par2name(yaml_data['multiphase']['phase0']['adimensional']['mu'])
    if(name_mult0_params != ''):
        name_mult0_params = name_mult0_params[1:]
    if(get_param_or_default(yaml_data, ['multiphase', 'phase1', 'adimensional', 'rho'], 1.0) != 1.0):
        name_mult1_params += '_\u03C1' + par2name(yaml_data['multiphase']['phase1']['adimensional']['rho'])
    if(get_param_or_default(yaml_data, ['multiphase', 'phase1', 'adimensional', 'mu'], 1.0) != 1.0):
        name_mult1_params += '_\u03BC' + par2name(yaml_data['multiphase']['phase1']['adimensional']['mu'])
    if(name_mult1_params != ''):
        name_mult1_params = name_mult1_params[1:]

    name_flowtype0 = ''
    name_flowtype0_model = ''
    name_flowtype0_params = ''
    name_flowtype1 = ''
    name_flowtype1_model = ''
    name_flowtype1_params = ''
    name_eoflow0 = ''
    name_eoflow0_model = ''
    name_eoflow0_params = ''
    name_eoflow1 = ''
    name_eoflow1_model = ''
    name_eoflow1_params = ''

    if(check_and_compare(yaml_data, ['multiphase','contr','eoflow0'], True, False)):
        yaml_eo = yaml_data['multiphase']['electroosmotic']
        name_eoflow0 = name_eo_symmetry+'EO'
        if yaml_eo['contr']['model'] == 'pb':
            name_eoflow0_model = 'pb'
        elif yaml_eo['contr']['model'] == 'pbdh':
            name_eoflow0_model = 'pbdh'
        else:
            name_eoflow0_model = ''  # pnp needs no description
        param_path = yaml_eo['phase0']['adimensional']
        name_eoflow0_params = '\u03BA' + par2name(math.sqrt(2.0*param_path['alpha']*param_path['delta']))
        name_eoflow0_params += '_\u03B1' + par2name(param_path['alpha'])
        if yaml_eo['contr']['model'] == 'pnp':
            name_eoflow0_params += '_Pe' + par2name(yaml_data['multiphase']['electroosmotic']['phase0']['model_pnp']['Pe'])
        name_eoflow0_params += '_Ex' + par2name(param_path['Ex'])
        if(get_param_or_default(param_path, ['perm'], 1.0) != 1.0):
            name_eoflow0_params += '_\u03B5\u2091' + par2name(param_path['perm'])
    elif(yaml_data['multiphase']['contr']['flowtype0'] == 'newtonian'):
        name_flowtype0 = 'Newt'

    if(yaml_data['multiphase']['contr']['flowtype0'] == 'viscoelastic'):
        yaml_visc = yaml_data['multiphase']['viscoelastic']
        name_flowtype0 = 'Visc'
        name_flowtype0_model = yaml_visc['contr']['model0']
        name_flowtype0_params = 'De' + par2name(yaml_visc['phase0']['adimensional']['De']) + '_\u03B2' + par2name(yaml_visc['phase0']['adimensional']['beta'])
        if(name_flowtype0_model == 'giesekus'):
            name_flowtype0_params += '_\u03B1' + par2name(yaml_visc['phase0']['model_giesekus']['alpha'])
        elif(name_flowtype0_model == 'lptt'):
            name_flowtype0_params += '_\u03B5' + par2name(yaml_visc['phase0']['model_lptt']['epsilon']) 
            if(yaml_visc['phase0']['model_lptt']['xi'] != 0.0):
                name_flowtype0_params += '_\u03BE' + par2name(yaml_visc['phase0']['model_lptt']['xi'])
        elif(name_flowtype0_model == 'gptt'):
            name_flowtype0_params += '_\u03B5' + par2name(yaml_visc['phase0']['model_gptt']['epsilon']) 
            if(yaml_visc['phase0']['model_gptt']['xi'] != 0.0):
                name_flowtype0_params += '_\u03BE' + par2name(yaml_visc['phase0']['model_gptt']['xi'])
            if(yaml_visc['phase0']['model_gptt']['alpha_gptt'] != 1.0):
                name_flowtype0_params += '_\u03B1_ml' + par2name(yaml_visc['phase0']['model_gptt']['alpha_gptt'])
            if(yaml_visc['phase0']['model_gptt']['beta_gptt'] != 1.0):
                name_flowtype0_params += '_\u03B2_ml' + par2name(yaml_visc['phase0']['model_gptt']['beta_gptt'])
        elif(name_flowtype0_model == 'fene_p'):
            name_flowtype0_params += '_L\u00B2' + par2name(yaml_visc['phase0']['model_fene_p']['L2'])
        elif(name_flowtype0_model == 'e_fene'):
            name_flowtype0_params += '_L\u00B2' + par2name(yaml_visc ['phase0']['model_e_fene']['L2'])
            name_flowtype0_params += '_\u03BB' + par2name(yaml_visc['phase0']['model_e_fene']['lambda']) + '_E' + str(yaml_visc['phase0']['model_e_fene']['E'])

    if(check_and_compare(yaml_data, ['multiphase','contr','eoflow1'], True, False)):
        yaml_eo = yaml_data['multiphase']['electroosmotic']
        if(name_eoflow0 == ''):
            name_eoflow1 = name_eo_symmetry+'EO'
        else:
            name_eoflow1 = 'EO'
        if yaml_eo['contr']['model'] == 'pb':
            name_eoflow1_model = 'pb'
        elif yaml_eo['contr']['model'] == 'pbdh':
            name_eoflow1_model = 'pbdh'
        else:
            name_eoflow1_model = ''  # pnp needs no description
        param_path = yaml_eo['phase1']['adimensional']
        name_eoflow1_params = '\u03BA' + par2name(math.sqrt(2.0*param_path['alpha']*param_path['delta']))
        name_eoflow1_params += '_\u03B1' + par2name(param_path['alpha'])
        if yaml_eo['contr']['model'] == 'pnp':
            name_eoflow1_params += '_Pe' + par2name(yaml_data['multiphase']['electroosmotic']['phase1']['model_pnp']['Pe'])
        name_eoflow1_params += '_Ex' + par2name(param_path['Ex'])
        if(get_param_or_default(param_path, ['perm'], 1.0) != 1.0):
            name_eoflow1_params += '_\u03B5\u2091' + par2name(param_path['perm'])
    elif(yaml_data['multiphase']['contr']['flowtype1'] == 'newtonian'):
        name_flowtype1 = 'Newt'

    if(yaml_data['multiphase']['contr']['flowtype1'] == 'viscoelastic'):
        yaml_visc = yaml_data['multiphase']['viscoelastic']
        name_flowtype1 = 'Visc'
        name_flowtype1_model = yaml_visc['contr']['model1']
        name_flowtype1_params = 'De' + par2name(yaml_visc['phase1']['adimensional']['De']) + '_\u03B2' + par2name(yaml_visc['phase1']['adimensional']['beta'])
        if(name_flowtype1_model == 'giesekus'):
            name_flowtype1_params += '_\u03B1' + par2name(yaml_visc['phase1']['model_giesekus']['alpha'])
        elif(name_flowtype1_model == 'lptt'):
            name_flowtype1_params += '_\u03B5' + par2name(yaml_visc['phase1']['model_lptt']['epsilon']) 
            if(yaml_visc['phase1']['model_lptt']['xi'] != 0.0):
                name_flowtype1_params += '_\u03BE' + par2name(yaml_visc['phase1']['model_lptt']['xi'])
        elif(name_flowtype1_model == 'gptt'):
            name_flowtype1_params += '_\u03B5' + par2name(yaml_visc['phase1']['model_gptt']['epsilon']) 
            if(yaml_visc['phase1']['model_gptt']['xi'] != 0.0):
                name_flowtype1_params += '_\u03BE' + par2name(yaml_visc['phase1']['model_gptt']['xi'])
            if(yaml_visc['phase1']['model_gptt']['alpha_gptt'] != 1.0):
                name_flowtype1_params += '_\u03B1_ml' + par2name(yaml_visc['phase1']['model_gptt']['alpha_gptt'])
            if(yaml_visc['phase1']['model_gptt']['beta_gptt'] != 1.0):
                name_flowtype1_params += '_\u03B2_ml' + par2name(yaml_visc['phase1']['model_gptt']['beta_gptt'])
        elif(name_flowtype1_model == 'fene_p'):
            name_flowtype1_params += '_L\u00B2' + par2name(yaml_visc['phase1']['model_fene_p']['L2'])
        elif(name_flowtype1_model == 'e_fene'):
            name_flowtype1_params += '_L\u00B2' + par2name(yaml_visc['phase1']['model_e_fene']['L2'])
            name_flowtype1_params += '_\u03BB' + par2name(yaml_visc['phase1']['model_e_fene']['lambda']) + '_E' + par2name(yaml_visc['phase1']['model_e_fene']['E'])

    if name_flowtype0 != '':
        name_flowtype0 = name_flowtype0 + '-'
    if name_flowtype0_model != '':
        name_flowtype0_model = name_flowtype0_model + '-'
    if name_eoflow0 != '':
        name_eoflow0 = name_eoflow0 + '-'
    if name_eoflow0_model != '':
        name_eoflow0_model = name_eoflow0_model + '-'
    if name_mult0_params != '':
        name_mult0_params = '-' + name_mult0_params
    if name_eoflow0_params != '':
        name_eoflow0_params = '-' + name_eoflow0_params
    if name_flowtype0_params != '':
        name_flowtype0_params = '-' + name_flowtype0_params

    if name_flowtype1 != '':
        name_flowtype1 = name_flowtype1 + '-'
    if name_flowtype1_model != '':
        name_flowtype1_model = name_flowtype1_model + '-'
    if name_eoflow1 != '':
        name_eoflow1 = name_eoflow1 + '-'
    if name_eoflow1_model != '':
        name_eoflow1_model = name_eoflow1_model + '-'
    if name_mult1_params != '':
        name_mult1_params = '-' + name_mult1_params
    if name_eoflow1_params != '':
        name_eoflow1_params = '-' + name_eoflow1_params
    if name_flowtype1_params != '':
        name_flowtype1_params = '-' + name_flowtype1_params

    params0 = name_mult0_params + name_eoflow0_params + name_flowtype0_params
    if params0 != '':
        params0 = '-par0' + params0
    params1 = name_mult1_params + name_eoflow1_params + name_flowtype1_params
    if params1 != '':
        params1 = '-par1' + params1 
        if 'Newt' not in name_flowtype0:
            params1 = '-' + params1

    name = f"{name_phase}"
    name += f"{name_eoflow0}{name_eoflow0_model}{name_flowtype0}{name_flowtype0_model}"
    if name[-1] == '-':
        name = name[:-1]
    name += f"__{name_eoflow1}{name_eoflow1_model}{name_flowtype1}{name_flowtype1_model}"
    name += f"{params0}"
    name += f"{params1}"
    name += f"{name_mult_params}{name_params}"

print(name)


