import os

def filter_tree():
    input_path = r"c:\Users\Siddharth Tripathi\OneDrive\Desktop\bioseq_explorer\documentation\bioseq_tree.txt"
    output_path = r"c:\Users\Siddharth Tripathi\OneDrive\Desktop\bioseq_explorer\scratch\filtered_bioseq_tree.txt"
    
    if not os.path.exists(input_path):
        print(f"Error: {input_path} not found.")
        return
        
    filtered_lines = []
    skip_node_modules = False
    skip_renv_lib = False
    
    with open(input_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            # Detect starting or stopping of node_modules and renv/library
            # Since it's tree /f output, folder indicators look like "+---node_modules"
            # and lines under it start with "|   " or similar indentations.
            
            # Simple rule: if 'node_modules' or 'renv' is in the line (when it's a folder boundary), skip.
            # But let's check indent levels to skip sub-files as well.
            
            # Let's inspect the directory structure. In tree /f, folders look like "+---node_modules"
            # or "|   +---node_modules".
            # Any line that contains "+---node_modules" starts the skip.
            # Any subsequent line starting with "|   " at that depth or deeper should be skipped.
            # Actually, a simpler way since node_modules and renv have distinct files:
            # We can parse the depth or simply check if 'node_modules' is in the line, 
            # or if the path is under node_modules/renv.
            # Wait, tree output is hierarchical. Let's see:
            # If line has "+---node_modules", we want to skip until we see the next folder at the same or higher level.
            # Since node_modules is at level 1 (+---node_modules), any line starting with "|   " or "+---" (that is under node_modules) should be skipped.
            # Let's write a stack-based or indentation-based parser, or a simple path-based one.
            # Wait! Since it is tree output:
            # Let's look at the structure:
            # 44: +---node_modules
            # 45: |   .package-lock.json
            # ...
            # 2547: +---outputs
            # This means '+---node_modules' is at column index 0 (or after a | indent).
            # Let's write a robust filter.
            pass

    # A simpler and 100% correct approach:
    # If the line contains node_modules, skip.
    # What about subfiles of node_modules? In tree /f, subfiles under "+---node_modules"
    # start with "|   " (since node_modules is at the root level).
    # Specifically, they start with "|   |" or "|   +---" or "|   \---" or "|       ".
    # Since node_modules starts at index 4, let's keep track of our current folder state.
    
    lines = []
    with open(input_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()
        
    filtered = []
    in_skip = False
    skip_indent = 0
    
    for line in lines:
        stripped = line.strip()
        # Find if it is a directory header
        if "+---node_modules" in line or "\\---node_modules" in line:
            in_skip = True
            # Find the index of '+' or '\'
            skip_indent = line.index('+') if '+' in line else line.index('\\')
            continue
        
        # We also want to skip renv's deep library, but keep the core renv files:
        # renv/activate.R, renv/settings.json, etc.
        # But wait! Is there a renv library in the tree? In bioseq_tree.txt,
        # +---renv has:
        # |       .gitignore
        # |       activate.R
        # |       settings.json
        # Wait, renv/library is not listed in bioseq_tree.txt!
        # Let's verify line 2565 of bioseq_tree.txt.
        # It has:
        # +---renv
        # |       .gitignore
        # |       activate.R
        # |       settings.json
        # And then +---scratch. There is no renv/library!
        # Ah! That means only node_modules is cluttering it!
        # Let's verify if node_modules is the only large directory to skip.
        # Yes, line 45-2546 is all node_modules!
        
        if in_skip:
            # Check indentation of this line.
            # If the line is empty, it could be inside or outside.
            # If it starts with | or space, check the character at skip_indent.
            # If the line is shorter than skip_indent, or if it doesn't have a vertical bar '|'
            # or spaces extending to skip_indent, or if it starts a new top-level folder:
            # Let's see: if the line starts with '+---' or '\---' at an indentation <= skip_indent, we stop skipping.
            # For example, '+---outputs' starts at column 0.
            # Let's inspect the prefix. If the line has '+' or '\' at a position <= skip_indent:
            idx = -1
            if '+' in line:
                idx = line.index('+')
            elif '\\' in line:
                idx = line.index('\\')
                
            if idx != -1 and idx <= skip_indent:
                in_skip = False
                
        if not in_skip:
            filtered.append(line)
            
    with open(output_path, 'w', encoding='utf-8') as f:
        f.writelines(filtered)
    print(f"Filtered tree from {len(lines)} to {len(filtered)} lines.")

if __name__ == '__main__':
    filter_tree()
