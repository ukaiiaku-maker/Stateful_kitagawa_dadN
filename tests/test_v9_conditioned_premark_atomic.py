import copy, math
import hashlib
import numpy as np

from arrhenius_fracture.v9_conditioned_premark import ConditionedBranchStreams, commit_conditioned_attempt


class DummyClock:
    def __init__(self):
        self.global_cumulative_action=math.log(2);self.global_threshold_action=math.log(2);self.attempt_count=0
        self._mark_rng=np.random.default_rng(1);self._hazard_rng=np.random.default_rng(2)
    def copy(self): return copy.deepcopy(self)
    def select_mark(self,nodes,available,logs,weights,**kw):
        ids,p,_=__import__('arrhenius_fracture.v9_pd_shared_root_marked_cleavage',fromlist=['normalized_available_site_marks']).normalized_available_site_marks(nodes,available,logs,weights)
        j=int(self._mark_rng.choice(len(ids),p=p)); site=int(ids[j]); self.attempt_count+=1
        self.global_threshold_action=self.global_cumulative_action+float(self._hazard_rng.exponential())
        return dict(site_id=site,pd_node_id=int(nodes[site]),mark_probability=float(p[j]),**kw)

class State:
    def __init__(self):
        self.site_status=np.zeros(3,np.uint8);self.site_node_index=np.array([0,1,1]);self.birth=[]
        self.candidate_sites=np.array([1,2])
        self.site_transition_threshold=np.full(3,99.);self.site_transition_cumulative_hazard=np.full(3,88.)
        self.site_transition_outcome_uniform=np.full(3,.5)
class Patch:
    initiation_weight=np.array([1.,2.])
    xy=np.array([[0.,0.],[1.,0.]])
    def create_marked_embryo(self,s,site,cycle): s.site_status[site]=1;s.birth.append((site,cycle))

def capsule():
    from arrhenius_fracture.v9_pd_shared_root_marked_cleavage import normalized_available_site_marks
    ids,p,_=normalized_available_site_marks([0,1,1],[1,1,1],[0.,.4],[1.,2.])
    candidates=np.array([1,2]); geometry=Patch.xy
    return dict(protocol='conditioned_attempt_premark',conditioned_first_attempt_pending=True,
      topology_continuation_permitted=True,exact_phase_localization_present=True,
      physical_threshold_draw=False,mark_rng_consumed=False,conditioned_action=math.log(2),
      candidate_population_sha256=hashlib.sha256(np.ascontiguousarray(candidates).view(np.uint8)).hexdigest(),
      geometry_sha256=hashlib.sha256(np.ascontiguousarray(geometry).view(np.uint8)).hexdigest(),
      branch_seed_namespace='test',source_analysis_checkpoint='a',source_replay_checkpoint='b',
      localization=dict(crossing_cycle=12.5,phase_index=2,phase_fraction=.25,
        local_cleavage_log_propensity_s=[0.,.4],local_opening_stress_Pa=[1.,2.],
        local_backstress_Pa=[3.,4.],local_state_shift_eV=[5.,6.],
        available_site_ids=ids.tolist(),normalized_mark_probability=p.tolist(),mark_entropy_nats=1.))

def test_atomic_conditioned_commit_is_restart_deterministic_and_does_not_mutate_source():
    c=DummyClock();s=State(); before=copy.deepcopy((c.__dict__,s.__dict__))
    a=commit_conditioned_attempt(c,Patch(),s,capsule(),'B003')
    b=commit_conditioned_attempt(c,Patch(),s,capsule(),'B003')
    assert a['event']==b['event'] and a['stream_seeds']==b['stream_seeds']
    assert a['clock'].attempt_count==1 and sum(a['pd_state'].site_status==1)==1
    site=a['event']['site_id']
    assert a['pd_state'].site_transition_cumulative_hazard[site]==0
    assert a['pd_state'].site_transition_threshold[site]!=99
    assert c.attempt_count==0 and not np.any(s.site_status) and before[0]['global_threshold_action']==c.global_threshold_action

def test_atomic_conditioned_commit_rolls_back_on_invalid_distribution():
    c=DummyClock();s=State();bad=capsule();bad['localization']['normalized_mark_probability'][0]+=.1
    try: commit_conditioned_attempt(c,Patch(),s,bad,'B000')
    except RuntimeError: pass
    else: raise AssertionError('invalid mark replay did not fail closed')
    assert c.attempt_count==0 and not np.any(s.site_status)

def test_three_stream_identities_are_independently_controllable():
    cap=capsule(); c=DummyClock(); s=State(); p=Patch()
    base=ConditionedBranchStreams('base','m','t','r')
    a=commit_conditioned_attempt(c,p,s,cap,base)
    t=commit_conditioned_attempt(c,p,s,cap,ConditionedBranchStreams('t','m','t2','r'))
    r=commit_conditioned_attempt(c,p,s,cap,ConditionedBranchStreams('r','m','t','r2'))
    assert a['event']['site_id']==t['event']['site_id']
    assert a['selected_site_transition_threshold']!=t['selected_site_transition_threshold']
    assert a['event']==r['event']
    assert a['selected_site_transition_threshold']==r['selected_site_transition_threshold']
    assert a['clock'].global_threshold_action!=r['clock'].global_threshold_action
