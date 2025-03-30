package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

// saveIndex writes the suffix entries to a file (one entry per line: pos, line, LCP).
func saveIndex(filename string, entries []SuffixEntry) error {
	var lines []string
	for _, entry := range entries {
		line := fmt.Sprintf("%d %d %d", entry.Pos, entry.Line, entry.LCP)
		lines = append(lines, line)
	}
	content := strings.Join(lines, "\n")
	return os.WriteFile(filename, []byte(content), 0644)
}

// loadIndex reads the suffix entries from a file.
func loadIndex(filename string) ([]SuffixEntry, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	var entries []SuffixEntry
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		parts := strings.Fields(line)
		if len(parts) != 3 {
			continue
		}
		pos, err1 := strconv.Atoi(parts[0])
		lineNum, err2 := strconv.Atoi(parts[1])
		lcpVal, err3 := strconv.Atoi(parts[2])
		if err1 != nil || err2 != nil || err3 != nil {
			continue
		}
		entries = append(entries, SuffixEntry{Pos: pos, Line: lineNum, LCP: lcpVal})
	}
	if err := scanner.Err(); err != nil {
		return nil, err
	}
	return entries, nil
}
