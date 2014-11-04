from AXLmodelSpat import *
import bayessb


if __name__ == '__main__':
    opts = bayessb.MCMCOpts()
    opts.model = Model();
    opts.estimate_params = opts.model.params()
    opts.likelihood_fn = likelihood
    opts.nsteps = int(1E6)
    opts.step_fn = step;
    opts.prior_fn = prior;
    opts.use_hessian = 1;
    opts.anneal_length = int(3E5);
    
    global mcmc
    
    mcmc = bayessb.MCMC(opts)
    mcmc.run()
    
