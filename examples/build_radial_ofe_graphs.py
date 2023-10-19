from openfe import SolventComponent
from openfe.setup import RHFEAlchemicalNetworkPlanner
from openfe.setup.ligand_network_planning import generate_radial_network


mapper = LomapAtomMapper()
network_planner = lambda ligands, mappers, scorer : generate_radial_network(ligands=ligands, 
                                                                            central_ligand=components[-5],
                                                                            mappers=mappers, scorer=scorer)

alchem_map_plan = RHFEAlchemicalNetworkPlanner(name="radial_RHFE_lomap", mappers=[mapper], mapping_scorer=None,
                             ligand_network_planner=network_planner)


solvent=SolventComponent()
alchemical_map_lomap = alchem_map_plan(ligands=rest_components, solvent=solvent)
