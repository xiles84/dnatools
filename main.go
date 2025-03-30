package main

import (
	"flag"
	"fmt"
	"os"
	"strings"
)

// SuffixEntry holds the suffix array entry, the originating line, and LCP value.
type SuffixEntry struct {
	Pos  int
	Line int
	LCP  int
}

func main() {
	runApp(os.Args[1:])
}

func runApp(args []string) {
	fs := flag.NewFlagSet("dnatools", flag.ExitOnError)
	indexMode := fs.Bool("m", false, "Index mode: build suffix array index")
	searchQueryStr := fs.String("s", "", "Search mode: search for sequence using suffix array")
	trieFile := fs.String("t", "", "Trie search mode: file containing multiple query patterns (one per line)")
	fileName := fs.String("f", "genoma.txt", "Genome file name")
	fs.Parse(args)

	// Read genome file.
	data, err := os.ReadFile(*fileName)
	if err != nil {
		fmt.Println("Error reading genome file:", err)
		os.Exit(1)
	}

	// Build genome from file lines. Each nonempty line is a DNA sequence.
	lines := strings.Split(string(data), "\n")
	var sequences []string
	for _, line := range lines {
		trimmed := strings.TrimSpace(line)
		if trimmed != "" {
			sequences = append(sequences, trimmed)
		}
	}
	var genomeBuilder strings.Builder
	var lineMap []int
	for i, seq := range sequences {
		for _, ch := range seq {
			genomeBuilder.WriteRune(ch)
			lineMap = append(lineMap, i) // record originating line index
		}
		// Append a separator if not the last sequence.
		if i < len(sequences)-1 {
			genomeBuilder.WriteByte('$')
			lineMap = append(lineMap, -1) // -1 indicates separator
		}
	}
	genome := genomeBuilder.String()

	// Index mode using suffix array (with LCP) remains as before.
	if *indexMode {
		fmt.Println("Building suffix array index using SAIS algorithm...")
		encoded, alphabetSize := encodeString(genome)
		sa := SAISEntryPoint(encoded, alphabetSize)
		lcp := computeLCP(genome, sa)
		entries := make([]SuffixEntry, len(sa))
		for i, pos := range sa {
			lineNum := -1
			if pos < len(lineMap) {
				lineNum = lineMap[pos]
			}
			entries[i] = SuffixEntry{Pos: pos, Line: lineNum, LCP: lcp[i]}
		}
		err = saveIndex("sa.idx", entries)
		if err != nil {
			fmt.Println("Error saving index:", err)
			os.Exit(1)
		}
		fmt.Println("Index built and saved to sa.idx")
		// Suffix array search mode.
	} else if *searchQueryStr != "" {
		entries, err := loadIndex("sa.idx")
		if err != nil {
			fmt.Println("Error loading index:", err)
			os.Exit(1)
		}
		fmt.Printf("Searching for sequence: %s\n", *searchQueryStr)
		results := searchSequence(genome, entries, *searchQueryStr)
		if len(results) == 0 {
			fmt.Println("Sequence not found.")
		} else {
			fmt.Println("Sequence found at positions (global position, DNA line):")
			for _, entry := range results {
				fmt.Printf("(%d, %d) ", entry.Pos, entry.Line)
			}
			fmt.Println()
		}
		// Trie search mode: used with the -t flag.
	} else if *trieFile != "" {
		// Read the file containing multiple query patterns.
		patternData, err := os.ReadFile(*trieFile)
		if err != nil {
			fmt.Println("Error reading trie file:", err)
			os.Exit(1)
		}
		patternLines := strings.Split(string(patternData), "\n")
		var patterns []string
		for _, l := range patternLines {
			trimmed := strings.TrimSpace(l)
			if trimmed != "" {
				patterns = append(patterns, trimmed)
			}
		}
		// Build the trie from the query patterns.
		trie := NewTrie()
		for _, pat := range patterns {
			trie.Insert(pat)
		}
		// Search the genome using the trie.
		results := searchTrie(genome, trie)
		// Print results, annotating each found position with its DNA line from lineMap.
		for pat, positions := range results {
			var annotated []string
			for _, pos := range positions {
				lineNum := -1
				if pos < len(lineMap) {
					lineNum = lineMap[pos]
				}
				annotated = append(annotated, fmt.Sprintf("(%d, line %d)", pos, lineNum))
			}
			fmt.Printf("Pattern %q found at positions: %v\n", pat, annotated)
		}
	} else {
		fmt.Println("Please provide -m to build index, -s <sequence> for suffix array search, or -t <file> for trie search (used with -f).")
	}
}
