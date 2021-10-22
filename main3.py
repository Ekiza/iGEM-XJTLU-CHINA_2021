import pinetree as pt
# import random
import time


def execute(output, infnum):
    sim = pt.Model(cell_volume=1.1e-15)
    t4genome = pt.Genome(name="editedT4genome", length=169713)
    t4genomer = pt.Genome(name="editedT4genomer", length=169713)

    # define polymerase
    sim.add_polymerase(name="rnapol", copy_number=2500, speed=120, footprint=5)
    sim.add_polymerase(name="E_gp33", copy_number=0, speed=120, footprint=5)
    sim.add_polymerase(name="E_gp33_gp45", copy_number=0, speed=120, footprint=5)

    # add ribosome
    sim.add_ribosome(copy_number=27000, speed=63, footprint=5)

    # define promoters and terminators
    t4genome.add_promoter(name="pholin", start=161011, stop=161024, interactions={"E_gp33": 46, "E_gp33_gp45": 1.5e4})
    t4genome.add_terminator(name="tholin", start=161734, stop=161773, efficiency={"E_gp33": 1.0, "E_gp33_gp45": 1.0})

    t4genome.add_promoter(name="pluxr", start=104766, stop=104835, interactions={"E_gp33": 46, "E_gp33_gp45": 1.5e4})
    t4genome.add_terminator(name="tluxr", start=107347, stop=107377, efficiency={"E_gp33": 1.0, "E_gp33_gp45": 1.0})

    t4genomer.add_promoter(name="pgp33", start=18177, stop=18194, interactions={"rnapol": 6413.43})
    t4genomer.add_terminator(name="tgp33", start=21979, stop=22007, efficiency={"rnapol": 1.0})

    t4genomer.add_promoter(name="pgp45", start=137044, stop=137084, interactions={"rnapol": 1430.37})
    t4genomer.add_terminator(name="tgp45", start=137802, stop=137815, efficiency={"rnapol": 1.0})

    t4genomer.add_promoter(name="pgp55", start=128420, stop=128489, interactions={"rnapol": 1015.18})
    t4genomer.add_terminator(name="tgp55", start=130983, stop=131018, efficiency={"rnapol": 1.0})

    t4genomer.add_promoter(name="pri", start=109952, stop=109993, interactions={"E_gp33": 46, "E_gp33_gp45": 1.5e4})
    t4genomer.add_terminator(name="tri", start=110901, stop=110929, efficiency={"rnapol": 1.0})

    # define genes
    t4genome.add_gene(name="LuxR", start=105261, stop=106070, rbs_start=105264, rbs_stop=105286, rbs_strength=833392)
    t4genome.add_gene(name="holin", start=161031, stop=161687, rbs_start=161025, rbs_stop=161030, rbs_strength=21055)
    t4genomer.add_gene(name="gp33", start=19373, stop=19711, rbs_start=19348, rbs_stop=19365, rbs_strength=1163895)
    t4genomer.add_gene(name="gp45", start=137111, stop=137797, rbs_start=137085, rbs_stop=137110, rbs_strength=103120)
    t4genomer.add_gene(name="gp55", start=129557, stop=130114, rbs_start=129535, rbs_stop=129550, rbs_strength=963655)
    t4genomer.add_gene(name="RI", start=110219, stop=110512, rbs_start=110207, rbs_stop=110218, rbs_strength=174019)

    # define species
    sim.add_species("E", 0)
    sim.add_species("2T-2RI", 0)
    sim.add_species("2T", 0)
    sim.add_species("2RI", 0)
    sim.add_species("DNA", infnum * 60000)
    sim.add_species("2T-2RI-DNA", 0)

    # define reactions
    sim.add_reaction(1e6, ["rnapol", "gp55"], ["E"])
    sim.add_reaction(1, ["E"], ["rnapol", "gp55"])
    sim.add_reaction(1e6, ["E", "gp33"], ["E_gp33"])
    sim.add_reaction(1, ["E_gp33"], ["E", "gp33"])
    sim.add_reaction(1e6, ["E_gp33", "gp45"], ["E_gp33_gp45"])
    sim.add_reaction(1, ["E_gp33_gp45"], ["E_gp33", "gp45"])
    sim.add_reaction(5e4, ["holin", "holin"], ["2T"])
    sim.add_reaction(1e5, ["RI", "RI"], ["2RI"])
    sim.add_reaction(1e5, ["2T", "2RI"], ["2T-2RI"])
    sim.add_reaction(10, ["2T-2RI"], ["2T"])
    sim.add_reaction(1e4, ["2T-2RI", "DNA"], ["2T-2RI-DNA"])  # more phage DNA can make 2T-2RI tetramer stable
    sim.add_reaction(10, ["2T"], ["holin", "holin"])
    sim.add_reaction(10, ["2RI"], ["RI", "RI"])

    # run the simulation
    sim.register_genome(t4genome)
    sim.register_genome(t4genomer)
    # seed = random.randint(0, 2147483647)
    # sim.seed(seed)
    # outputstr = str(infnum) + "." + str(cellnum) + "." + str(p1) + "." + str(c1) + ".test.tsv"
    sim.seed(34)
    sim.simulate(time_limit=9000, time_step=500, output=output + str(infnum) + ".tsv")


if __name__ == "__main__":
    time_start = time.time()
    for num in range(16, 21):
        execute("real.", num)
        print(num)
    time_end = time.time()
    print('totally cost', time_end - time_start)
