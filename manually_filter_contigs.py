import pandas as pd
import tkinter as tk
from tkinter import messagebox
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
import argparse

# Load data from CSV files
def load_data(records_file, kmer_file, selection_file):
    #tabulated files 
    records_names =['Assembly', 'contig', 'haplotype', 'values']
    records = pd.read_csv(records_file, sep='\t', names=records_names, header=None)
    #print(records)
    
    kmer_names =['Assembly', 'pair', 'contig', 'values']
    kmer = pd.read_csv(kmer_file, sep='\t', names=kmer_names, header=None)  
    #print(kmer)

    #uniq file 
    selection_names = ['contig_name1','contig_name2']
    selection = pd.read_csv(selection_file, sep="\t", names=selection_names, header=None)
    #print(selection)
    #selection = False
    return records, kmer, selection

# Function to fetch joined contig names from the selection DataFrame
def get_joined_contig_names(selection_df):
    contig_names = selection_df.apply(lambda row: f"{row['contig_name1']}-{row['contig_name2']}", axis=1)
    return contig_names

# Function to filter records DataFrame by contig name
def get_records_by_name(records_df, contig_name):
    return records_df[records_df['contig'] == contig_name]

# Function to filter kmer DataFrame by contig name
def get_kmer_by_name(kmer_df, contig_name, selected_contig):
    return kmer_df[(kmer_df['contig'] == contig_name) & (kmer_df['pair'] == selected_contig)]

# Function to update and display the stacked bar plot for selected contig names
def on_contig_select(event):
    try:
        # Get the selected contig from the Listbox
        selected_index = listbox.curselection()[0]
        selected_contig = listbox.get(selected_index)

        # Split the selected contig into contig_name1 and contig_name2
        contig_name1, contig_name2 = selected_contig.split('-')

        # Filter records and kmer DataFrames for both contigs
        data1_records = get_records_by_name(records_df, contig_name1)
        data2_records = get_records_by_name(records_df, contig_name2)
        #print(data1_records,data2_records)
        #print("==>", contig_name1, contig_name2, selected_contig)
        data1_kmer = get_kmer_by_name(kmer_df, contig_name1, selected_contig)
        data2_kmer = get_kmer_by_name(kmer_df, contig_name2, selected_contig)
        #print("datakmer =>", data1_kmer, data2_kmer)

        # Clear the previous plots
        for ax in axs.flatten():
            ax.clear()

        # Prepare plot data for a stacked bar plot
        def prepare_plot_data(data):
            #print(data)
            values = {}
            l = 0
            
            for d in data[['haplotype','values']].iterrows() :
                values[d[1][0]] = d[1][1].split()
                l = len(d[1][1].split(","))
            
            df_values = pd.DataFrame.from_dict(values)
            locations = [ i for i in range(l)]
            categories = sorted(data['haplotype'].unique())           
            stacked_data = {cat: [] for cat in categories}

            for loc in locations:
                count = 0
                for cat in categories:
                    count = int(float(df_values[cat][0].split(",")[loc]))
                    stacked_data[cat].append(count)

            return locations, stacked_data
            
        # Prepare plot data for both contigs from records table
        unique_locations1, stacked_counts1 = prepare_plot_data(data1_records)
        #print("data1 done")
        unique_locations2, stacked_counts2 = prepare_plot_data(data2_records)
        #print("stacked_data : ",stacked_counts1, stacked_counts2)
        
        # Create a stacked bar plot for contig_name1 from records
        bottom1 = np.zeros(len(unique_locations1))
        #print("bottom1 ",bottom1)
        for cat, cat_counts in stacked_counts1.items():
            print("cat = ",unique_locations1, cat, cat_counts)
            axs[0, 0].bar(unique_locations1, cat_counts, width=1000, bottom=bottom1, label=str(cat))
            #print("bottom1 ",bottom1)
            bottom1 += np.array(cat_counts)

        axs[0, 0].set_title(f'Total : {contig_name1}', fontsize=12)
        axs[0, 0].set_xlabel('Location', fontsize=10)
        axs[0, 0].set_ylabel('Count', fontsize=10)

        # Create a stacked bar plot for contig_name2 from records
        bottom2 = np.zeros(len(unique_locations2))
        for cat, cat_counts in stacked_counts2.items():
            print("cat = ",unique_locations2, cat, cat_counts)
            axs[0, 1].bar(unique_locations2, cat_counts, width=1000, bottom=bottom2, label=str(cat))
            bottom2 += np.array(cat_counts)

        axs[0, 1].set_title(f'Total : {contig_name2}', fontsize=12)
        axs[0, 1].set_xlabel('Location', fontsize=10)
        axs[0, 1].set_ylabel('Count', fontsize=10)

        # Plotting for kmer data
        l = 0
        l1 = []
            
        print("data1_kmer", data1_kmer)    
        for d in data1_kmer[['contig','values']].iterrows() :
            l = len(d[1][1].split(","))
            l1 = d[1][1].split(",")
            
        #l = len(data1_kmer[['contig','values']].head(1)['values'][:1].split(","))
        locations1 = [ i for i in range(l)]
        counts1 = [int(float(i)) for i in l1]
        print("locations and counts 1 : ", locations1, counts1)
        #for i in locations1 :
        #    print(i, locations1[i], counts1[i])
            
        axs[1, 0].bar(locations1, counts1, width=1000)
        axs[1, 0].set_title(f'Kmer: {contig_name1}', fontsize=12)
        axs[1, 0].set_xlabel('Location', fontsize=10)
        axs[1, 0].set_ylabel('Count', fontsize=10)

        l = 0
        l1 = []
            
        for d in data2_kmer[['contig','values']].iterrows() :
            l = len(d[1][1].split(","))
            l1 = d[1][1].split(",")
            
        #l = len(data2_kmer[['contig','values']].head(1)['values'][0].split(","))
        locations2 = [ i for i in range(l)]
        counts2 = [int(float(i)) for i in l1]
        #for i in locations2 :
        #    print(i, locations2[i], counts2[i])
            
        axs[1, 1].bar(locations2, counts2, width=1000)
        axs[1, 1].set_title(f'Kmer: {contig_name2}', fontsize=12)
        axs[1, 1].set_xlabel('Location', fontsize=10)
        axs[1, 1].set_ylabel('Count', fontsize=10)

        # Synchronize the x-axis limits for both graphs
        all_locations = np.concatenate([unique_locations1, unique_locations2, locations1, locations2])
        #all_locations = np.concatenate([unique_locations1, unique_locations2])
        max_x = max(all_locations) if len(all_locations) > 0 else 1
        for ax in axs.flatten():
            ax.set_xlim(0, max_x)

        # Redraw the canvas with the updated plots
        axs[1, 1].legend()
        canvas.draw()

    except IndexError:
        # No selection made
        return

# Function to add or remove a red vertical line at the clicked position
def on_click(event):
    global vertical_line  # Use a global variable to track the vertical line
    x_coord = event.xdata

    if x_coord is not None:
        if vertical_line is not None:
            # If a vertical line already exists, remove it
            vertical_line.remove()
            vertical_line = None
            left_button.config(state=tk.DISABLED)
            right_button.config(state=tk.DISABLED)
        else:
            # Draw a new vertical line at the clicked position
            vertical_line = axs[0, 0].axvline(x=x_coord, color='red', linestyle='--')
            left_button.config(state=tk.NORMAL)
            right_button.config(state=tk.NORMAL)

        # Redraw the canvas to display the changes
        canvas.draw()

        # Print the selected contig and x-coordinate of the clicked position
        if listbox.curselection():
            selected_index = listbox.curselection()[0]
            selected_contig = listbox.get(selected_index)
            print(f"Contig: {selected_contig}, Clicked at x-coordinate: {x_coord:.2f}")

# Function to handle the left button click
def on_left_click():
    print("left")

# Function to handle the right button click
def on_right_click():
    print("right")

# Main Application
if __name__ == "__main__":
    # Argument parsing
    parser = argparse.ArgumentParser(description="Plot contig and kmer data.")
    parser.add_argument("--records", required=True, help="Path to the records CSV file.")
    parser.add_argument("--kmer", required=True, help="Path to the kmer CSV file.")
    parser.add_argument("--selection", required=True, help="Path to the selection CSV file.")
    
    args = parser.parse_args()

    # Load the data from the CSV files
    records_df, kmer_df, selection_df = load_data(args.records, args.kmer, args.selection)

    # Create Tkinter window
    root = tk.Tk()
    root.title("Contig Selection and Plotting")

    # Create Listbox for contig selection with scrollbar
    frame = tk.Frame(root)
    frame.pack(pady=10)

    scrollbar = tk.Scrollbar(frame)
    scrollbar.pack(side=tk.RIGHT, fill=tk.Y)

    listbox = tk.Listbox(frame, width=50, height=10, yscrollcommand=scrollbar.set)
    listbox.pack()

    scrollbar.config(command=listbox.yview)

    # Fetch and display contig names in the Listbox
    contig_names = get_joined_contig_names(selection_df)
    for name in contig_names:
        listbox.insert(tk.END, name)

    # Create matplotlib figure and axes
    fig, axs = plt.subplots(2, 2, figsize=(12, 8))
    vertical_line = None  # Track the vertical line

    # Create a canvas for the figure
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.get_tk_widget().pack()

    # Create buttons for left and right clicks
    left_button = tk.Button(root, text="Left Click", command=on_left_click, state=tk.DISABLED)
    left_button.pack(side=tk.LEFT, padx=(20, 10))

    right_button = tk.Button(root, text="Right Click", command=on_right_click, state=tk.DISABLED)
    right_button.pack(side=tk.LEFT)

    # Bind selection and click events
    listbox.bind("<<ListboxSelect>>", on_contig_select)
    canvas.mpl_connect("button_press_event", on_click)

    # Start the Tkinter main loop
    root.mainloop()

