# Title anno_manhattan_plot_marker_per_probe

Title anno_manhattan_plot_marker_per_probe

## Usage

``` r
anno_manhattan_plot_marker_per_probe(
  probe_name_max = "cg11680158",
  probe_name_min = "cg11680158",
  max_sample = 0,
  min_sample = 0,
  min_signal_probe = 0,
  label_font_size = 3,
  hyper_color = "blue",
  hypo_color = "orange",
  non_outlier_color = "grey",
  limit_label_color = "blue",
  limit_line_color = "red",
  limit_line_color_median = "black",
  reference_samples_color = "cyan",
  show_labels = FALSE,
  result_folder = result_folder,
  maxResources = maxResources,
  parallel_strategy = parallel_strategy,
  ...
)
```

## Arguments

- probe_name_max:

  cg name of the probe to represent tyh probe with maximum burden

- probe_name_min:

  cg name of the probe to represent tyh probe with miniumal burden

- max_sample:

  max number of samples to plot

- min_sample:

  min number of samples to plot

- min_signal_probe:

  min signal value of the probe to be plotted

- label_font_size:

  size of the labels

- hyper_color:

  color to assign to hypermethylated probes

- hypo_color:

  color to assign to hypomethylated probes

- non_outlier_color:

  color to assign to probes that are not outliers

- limit_label_color:

  color to assign to the labels of the limit lines

- limit_line_color:

  color to assign to the limit lines

- limit_line_color_median:

  color to assign to the median limit line

- reference_samples_color:

  color to assign to the reference samples

- show_labels:

  show labels in the plot

- result_folder:

  foder where the results are stored

- maxResources:

  percentage of max system's resource to use

- parallel_strategy:

  strategy to use for parallelization

- ...:

  other parameter

## Value

Invisibly `NULL`. The function saves PNG plot files to disk as a side
effect.
