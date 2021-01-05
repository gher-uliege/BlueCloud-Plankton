The ipynb files are generated from julia by using:

```
jupytext  --to notebook DIVAndNN_analysis.jl
jupytext  --to notebook DIVAndNN_plot_res.jl
```

Do not edit the ipynb files without updating the scripts too.

```
jupytext --sync DIVAndNN_analysis.ipynb
jupytext --sync DIVAndNN_plot_res.ipynb
```
