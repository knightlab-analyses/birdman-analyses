import seaborn as sns

tool_list = ["birdman", "deseq2", "aldex2", "ancombc"]
tool_text_styling_dict = {
    "birdman": "BIRDMAn",
    "deseq2": "DESeq2",
    "aldex2": "ALDEx2",
    "ancombc": "ANCOM-BC",
}

tool_palette = dict(zip(tool_text_styling_dict.values(), sns.color_palette("Dark2", 5)))