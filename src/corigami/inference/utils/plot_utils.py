import os
import numpy as np
import pandas as pd

class MatrixPlot:

    def __init__(self, output_path, image, prefix, celltype, chr_name, start_pos):
        self.output_path = output_path,
        self.prefix = prefix
        self.celltype = celltype
        self.chr_name = chr_name
        self.start_pos = start_pos

        self.create_save_path(output_path, celltype, prefix)
        self.image = self.preprocess_image(image)

    def get_colormap(self):
        from matplotlib.colors import LinearSegmentedColormap
        color_map = LinearSegmentedColormap.from_list("bright_red", [(1,1,1),(1,0,0)])
        return color_map

    def create_save_path(self, output_path, celltype, prefix):
        self.save_path = f'{output_path}/{celltype}/{prefix}'
        os.makedirs(f'{self.save_path}/imgs', exist_ok = True)
        os.makedirs(f'{self.save_path}/npy', exist_ok = True)

    def preprocess_image(self, image):
        return image

    def plot(self, vmin = 0, vmax = 5):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize = (5, 5))
        color_map = self.get_colormap()
        ax.imshow(self.image, cmap = color_map, vmin = vmin, vmax = vmax)
        self.reformat_ticks(plt)
        return 

    def reformat_ticks(self, plt):
        # Rescale tick labels
        current_ticks = np.arange(0, 250, 50) / 0.8192
        plt.xticks(current_ticks, self.rescale_coordinates(current_ticks, self.start_pos))
        plt.yticks(current_ticks, self.rescale_coordinates(current_ticks, self.start_pos))
        # Format labels
        plt.ylabel('Genomic position (Mb)')
        plt.xlabel(f'Chr{self.chr_name.replace("chr", "")}: {self.start_pos} - {self.start_pos + 2097152} ')
        self.save_data(plt)

    def rescale_coordinates(self, coords, zero_position):
        scaling_ratio = 8192
        replaced_coords = coords * scaling_ratio + zero_position
        coords_mb = replaced_coords / 1000000
        str_list = [f'{item:.2f}' for item in coords_mb]
        return str_list

    def save_data(self, plt):
        plt.savefig(f'{self.save_path}/imgs/{self.chr_name}_{self.start_pos}.png', bbox_inches = 'tight')
        plt.close()
        np.save(f'{self.save_path}/npy/{self.chr_name}_{self.start_pos}', self.image)

class MatrixPlotDeletion(MatrixPlot):
    def __init__(self, output_path, image, prefix, celltype, chr_name, start_pos, deletion_start, deletion_width, padding_type, show_deletion_line = False):
        super().__init__(output_path, image, prefix, celltype, chr_name, start_pos)
        self.deletion_start = deletion_start
        self.deletion_width = deletion_width
        self.show_deletion_line = show_deletion_line
        self.padding_type = padding_type

    def reformat_ticks(self, plt):
        # Rescale tick labels
        breakpoint_start = (self.deletion_start - self.start_pos) / 10000 
        breakpoint_end = (self.deletion_start - self.start_pos + self.deletion_width) / 10000 
        # Used for generating ticks until the end of the window
        total_window_size = (self.deletion_width + 2097152 ) / 10000
        # Generate ticks before and after breakpoint
        before_ticks = np.arange(0, breakpoint_start - 50, 50) / 0.8192
        after_ticks = (np.arange((breakpoint_end // 50 + 2) * 50, total_window_size, 50) - self.deletion_width / 10000) / 0.8192
        breakpoint_locus = breakpoint_start / 0.8192
        # Actual coordinates for each tick
        current_ticks = np.append(before_ticks, after_ticks)
        current_ticks = np.append(current_ticks, breakpoint_start / 0.8192)
        # Genomic coordinates used for display location after deletion
        display_ticks = np.append(before_ticks, after_ticks + self.deletion_width / 10000 / 0.8192)
        display_ticks = np.append(display_ticks, breakpoint_start / 0.8192)
        if self.show_deletion_line:
            plt.axline((breakpoint_locus, 0), (breakpoint_locus, 209), c = 'black', alpha = 0.5)
            plt.axline((0, breakpoint_locus), (209, breakpoint_locus), c = 'black', alpha = 0.5)
        # Generate tick label text
        ticks_label = self.rescale_coordinates(display_ticks, self.start_pos)
        plt.yticks(current_ticks, ticks_label)
        ticks_label[-1] = f"{(self.deletion_start / 1000000):.2f}({(self.deletion_start + self.deletion_width) / 1000000:.2f})"
        plt.xticks(current_ticks, ticks_label)
        # Format labels
        plt.ylabel('Genomic position (Mb)')
        end_pos = self.start_pos + 2097152 + self.deletion_width
        plt.xlabel(f'Chr{self.chr_name.replace("chr", "")}: {self.start_pos} - {self.deletion_start} and {self.deletion_start + self.deletion_width} - {end_pos} ')
        self.save_data(plt)

    def save_data(self, plt):
        plt.savefig(f'{self.save_path}/imgs/{self.chr_name}_{self.start_pos}_del_{self.deletion_start}_{self.deletion_width}_padding_{self.padding_type}.png', bbox_inches = 'tight')
        plt.close()
        np.save(f'{self.save_path}/npy/{self.chr_name}_{self.start_pos}_del_{self.deletion_start}_{self.deletion_width}_padding_{self.padding_type}', self.image)

class MatrixPlotPointScreen(MatrixPlotDeletion):

    def plot(self, vmin = -1, vmax = 1):
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize = (5, 5))
        ax.imshow(self.image, cmap = 'RdBu_r', vmin = vmin, vmax = vmax)
        self.reformat_ticks(plt)
        return 

    def save_data(self, plt):
        plt.savefig(f'{self.save_path}/imgs/{self.chr_name}_{self.start_pos}_del_{self.deletion_start}_{self.deletion_width}_padding_{self.padding_type}_diff.png', bbox_inches = 'tight')
        plt.close()
        np.save(f'{self.save_path}/npy/{self.chr_name}_{self.start_pos}_del_{self.deletion_start}_{self.deletion_width}_padding_{self.padding_type}_diff', self.image)

class MatrixPlotScreen(MatrixPlot):
    def __init__(self, output_path, perturb_starts, perturb_ends, impact_score, tensor_diff, tensor_pred, tensor_deletion, prefix, celltype, chr_name, screen_start, screen_end, perturb_width, step_size, plot_impact_score):
        super().__init__(output_path, impact_score, prefix, celltype, chr_name, start_pos = None)
        self.perturb_starts = perturb_starts
        self.perturb_ends = perturb_ends
        self.impact_score = impact_score
        self.tensor_diff = tensor_diff
        self.tensor_pred = tensor_pred
        self.tensor_deletion = tensor_deletion
        self.screen_start = screen_start
        self.screen_end = screen_end
        self.perturb_width = perturb_width
        self.step_size = step_size
        self.plot_impact_score = plot_impact_score

    def create_save_path(self, output_path, celltype, prefix):
        self.save_path = f'{output_path}/{celltype}/{prefix}'
        os.makedirs(f'{self.save_path}/imgs', exist_ok = True)
        os.makedirs(f'{self.save_path}/npy', exist_ok = True)
        os.makedirs(f'{self.save_path}/bedgraph', exist_ok = True)

    def plot(self, vmin = -1, vmax = 1):
        import matplotlib.pyplot as plt
        height = 3
        width = 1 * np.log2(len(self.impact_score))
        fig, ax = plt.subplots(figsize = (width, height))
        self.plot_track(ax, self.impact_score, self.screen_start, self.step_size)
        self.reformat_ticks(plt)
        return plt

    def reformat_ticks(self, plt):
        # Format labels
        plt.xlabel('Genomic position (Mb)')

    def save_data(self, plt, save_pred, save_deletion, save_diff, save_impact_score, save_bedgraph):
        if self.plot_impact_score:
            plt.savefig(f'{self.save_path}/imgs/{self.chr_name}_screen_{self.screen_start}_{self.screen_end}_width_{self.perturb_width}_step_{self.step_size}.png', bbox_inches = 'tight')
            plt.close()
        if save_pred:
            np.save(f'{self.save_path}/npy/{self.chr_name}_screen_{self.screen_start}_{self.screen_end}_width_{self.perturb_width}_step_{self.step_size}_pred', self.tensor_pred)
        if save_deletion:
            np.save(f'{self.save_path}/npy/{self.chr_name}_screen_{self.screen_start}_{self.screen_end}_width_{self.perturb_width}_step_{self.step_size}_perturbed', self.tensor_deletion)
        if save_diff:
            np.save(f'{self.save_path}/npy/{self.chr_name}_screen_{self.screen_start}_{self.screen_end}_width_{self.perturb_width}_step_{self.step_size}_diff', self.tensor_diff)
        if save_impact_score:
            np.save(f'{self.save_path}/npy/{self.chr_name}_screen_{self.screen_start}_{self.screen_end}_width_{self.perturb_width}_step_{self.step_size}_impact_score', self.impact_score)
        if save_bedgraph:
            bedgraph_path = f'{self.save_path}/bedgraph/{self.chr_name}_screen_{self.screen_start}_{self.screen_end}_width_{self.perturb_width}_step_{self.step_size}_impact_score.bedgraph'
            self.save_bedgraph(self.chr_name, self.perturb_starts, self.perturb_ends, self.impact_score, bedgraph_path)

    def plot_track(self, ax, data, start, step):
        x = (np.array(range(len(data))) + int(start / step)) * step / 1000000
        width = min(self.perturb_width, int(step * 0.9)) / 1000000
        ax.bar(x, data, width = width)
        ax.margins(x=0)
        #ax.get_xaxis().set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        #ax.spines['bottom'].set_visible(False)
        ax.set_ylabel('Impact score')
        #ax.set_ylim(-1, 8)

    def save_bedgraph(self, chr_name, starts, ends, scores, output_file):
        df = pd.DataFrame({'chr': chr_name, 'start': starts, 'end': ends, 'score': scores})
        df.to_csv(output_file, sep = '\t', index = False, header = False)

