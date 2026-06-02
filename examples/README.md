# How to generate the benzene data

This is how I generated this file, I checked out this commit from openfe-benchmarks [`d5a027e4e3cb53e47d4b230e8ddffda274b70aff`](https://github.com/OpenFreeEnergy/openfe-benchmarks/tree/d5a027e4e3cb53e47d4b230e8ddffda274b70aff/openfe_benchmarks) and then from the `openfe_benchmarks` folder ran these commands:

```
>>> import benzenes
>>> components = benzenes.get_system().ligand_components
>>> not_lig = ["lig_4", "lig_7", "lig_2", "lig_3"]
>>> components = [c for c in components if (c.name not in not_lig)]
>>> for comp in components:
...     comp.to_msgpack(f"{comp.name}.msgpack")
...
>>>
```

I then tar'ed all the msgpack files into a gziped tarball. 
