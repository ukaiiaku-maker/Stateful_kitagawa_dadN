import numpy as np
from arrhenius_fracture.v9_pd_shared_root_marked_cleavage import normalized_available_site_marks

def test_conditioned_premark_probabilities_are_normalized_and_colocated_sites_split_equally():
    nodes=np.array([2,2,3]);available=np.ones(3,dtype=bool);logp=np.array([0.,0.,1.,0.])
    ids,p,_=normalized_available_site_marks(nodes,available,logp,np.ones(4))
    assert ids.tolist()==[0,1,2]
    assert np.isclose(p.sum(),1.)
    assert p[0]==p[1]
