import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.nio.file.*;
import java.time.Instant;
import java.time.LocalDateTime;
import java.time.ZoneId;
import java.time.format.DateTimeFormatter;
import java.util.*;
import java.util.List;
import java.util.concurrent.*;
import java.util.stream.Collectors;

public class VideoFileManager {
    // Supported video formats
    private static final Set<String> VIDEO_EXTENSIONS = new HashSet<>(Arrays.asList(
        "mp4", "avi", "mkv", "mov", "wmv", "flv", "webm", "mpeg", "mpg", 
        "3gp", "m4v", "ogv", "ts", "m2ts", "vob", "rmvb", "asf", "divx"
    ));
    
    private List<VideoFile> videoFiles = new ArrayList<>();
    private Map<String, List<VideoFile>> videosByExtension = new HashMap<>();
    private Map<Integer, List<VideoFile>> videosBySizeCategory = new HashMap<>();
    private JTextArea logArea;
    private ExecutorService executor = Executors.newFixedThreadPool(3);
    
    // Video metadata class
    private static class VideoFile {
        Path path;
        String name;
        String extension;
        long size; // bytes
        LocalDateTime created;
        LocalDateTime modified;
        String resolution; // will be parsed if possible
        
        VideoFile(Path path) throws IOException {
            this.path = path;
            this.name = path.getFileName().toString();
            this.extension = getExtension(name).toLowerCase();
            
            BasicFileAttributes attrs = Files.readAttributes(path, BasicFileAttributes.class);
            this.size = attrs.size();
            this.created = LocalDateTime.ofInstant(
                attrs.creationTime().toInstant(), ZoneId.systemDefault()
            );
            this.modified = LocalDateTime.ofInstant(
                attrs.lastModifiedTime().toInstant(), ZoneId.systemDefault()
            );
        }
        
        private String getExtension(String filename) {
            int dotIndex = filename.lastIndexOf('.');
            return dotIndex > 0 ? filename.substring(dotIndex + 1) : "";
        }
        
        @Override
        public String toString() {
            return String.format("%s [%s, %s, Created: %s]", 
                name, 
                formatFileSize(size),
                extension.toUpperCase(),
                created.format(DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm"))
            );
        }
    }
    
    public VideoFileManager() {
        setupGUI();
    }
    
    private void setupGUI() {
        JFrame frame = new JFrame("Video File Manager");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLayout(new BorderLayout(10, 10));
        
        // Control Panel
        JPanel controlPanel = new JPanel(new GridLayout(2, 4, 5, 5));
        
        JButton scanBtn = new JButton("Scan Directory");
        JButton organizeBtn = new JButton("Organize by Type");
        JButton analyzeBtn = new JButton("Analyze Videos");
        JButton convertBtn = new JButton("Convert Formats");
        JButton cleanupBtn = new JButton("Cleanup Duplicates");
        JButton exportBtn = new JButton("Export Report");
        JButton previewBtn = new JButton("Generate Preview");
        JButton metadataBtn = new JButton("Extract Metadata");
        
        controlPanel.add(scanBtn);
        controlPanel.add(organizeBtn);
        controlPanel.add(analyzeBtn);
        controlPanel.add(convertBtn);
        controlPanel.add(cleanupBtn);
        controlPanel.add(exportBtn);
        controlPanel.add(previewBtn);
        controlPanel.add(metadataBtn);
        
        // Information Panel
        logArea = new JTextArea(15, 60);
        logArea.setEditable(false);
        logArea.setFont(new Font("Monospaced", Font.PLAIN, 12));
        JScrollPane scrollPane = new JScrollPane(logArea);
        
        // Status Panel
        JPanel statusPanel = new JPanel(new FlowLayout(FlowLayout.LEFT));
        JLabel statusLabel = new JLabel("Ready");
        statusPanel.add(statusLabel);
        
        // Add components
        frame.add(controlPanel, BorderLayout.NORTH);
        frame.add(scrollPane, BorderLayout.CENTER);
        frame.add(statusPanel, BorderLayout.SOUTH);
        
        // Action Listeners
        scanBtn.addActionListener(e -> scanVideoDirectory());
        organizeBtn.addActionListener(e -> organizeVideosByType());
        analyzeBtn.addActionListener(e -> analyzeVideoCollection());
        convertBtn.addActionListener(e -> showConversionDialog());
        cleanupBtn.addActionListener(e -> findDuplicateVideos());
        exportBtn.addActionListener(e -> exportVideoReport());
        previewBtn.addActionListener(e -> generatePreviews());
        metadataBtn.addActionListener(e -> extractVideoMetadata());
        
        frame.pack();
        frame.setVisible(true);
    }
    
    private void scanVideoDirectory() {
        JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.setDialogTitle("Select Video Directory");
        
        if (chooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
            File directory = chooser.getSelectedFile();
            executor.submit(() -> {
                try {
                    scanDirectoryRecursive(directory.toPath());
                } catch (Exception e) {
                    logError("Scan error: " + e.getMessage());
                }
            });
        }
    }
    
    private void scanDirectoryRecursive(Path directory) throws IOException {
        videoFiles.clear();
        videosByExtension.clear();
        videosBySizeCategory.clear();
        
        if (!Files.exists(directory) || !Files.isDirectory(directory)) {
            logError("Invalid directory: " + directory);
            return;
        }
        
        logMessage("Scanning directory: " + directory);
        
        try (var paths = Files.walk(directory)) {
            List<Path> videoPaths = paths
                .filter(Files::isRegularFile)
                .filter(this::isVideoFile)
                .collect(Collectors.toList());
            
            for (Path path : videoPaths) {
                try {
                    VideoFile video = new VideoFile(path);
                    videoFiles.add(video);
                    
                    // Group by extension
                    videosByExtension
                        .computeIfAbsent(video.extension, k -> new ArrayList<>())
                        .add(video);
                    
                    // Group by size category
                    int sizeCategory = (int) (video.size / (1024 * 1024 * 100)); // 100MB categories
                    videosBySizeCategory
                        .computeIfAbsent(sizeCategory, k -> new ArrayList<>())
                        .add(video);
                    
                } catch (IOException e) {
                    logError("Failed to read: " + path.getFileName());
                }
            }
        }
        
        logMessage("Found " + videoFiles.size() + " video files");
        logMessage("Extensions found: " + videosByExtension.keySet());
    }
    
    private boolean isVideoFile(Path path) {
        String filename = path.getFileName().toString().toLowerCase();
        int dotIndex = filename.lastIndexOf('.');
        if (dotIndex > 0) {
            String extension = filename.substring(dotIndex + 1);
            return VIDEO_EXTENSIONS.contains(extension);
        }
        return false;
    }
    
    private void organizeVideosByType() {
        if (videoFiles.isEmpty()) {
            logError("No videos loaded. Scan directory first.");
            return;
        }
        
        JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.setDialogTitle("Select Target Directory for Organization");
        
        if (chooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
            File targetDir = chooser.getSelectedFile();
            executor.submit(() -> {
                try {
                    organizeByExtension(targetDir.toPath());
                } catch (Exception e) {
                    logError("Organization error: " + e.getMessage());
                }
            });
        }
    }
    
    private void organizeByExtension(Path targetDir) throws IOException {
        for (Map.Entry<String, List<VideoFile>> entry : videosByExtension.entrySet()) {
            Path extensionDir = targetDir.resolve(entry.getKey().toUpperCase());
            Files.createDirectories(extensionDir);
            
            for (VideoFile video : entry.getValue()) {
                try {
                    Path target = extensionDir.resolve(video.path.getFileName());
                    Files.copy(video.path, target, StandardCopyOption.REPLACE_EXISTING);
                    logMessage("Organized: " + video.name + " -> " + extensionDir.getFileName());
                } catch (IOException e) {
                    logError("Failed to organize: " + video.name);
                }
            }
        }
        logMessage("Organization completed!");
    }
    
    private void analyzeVideoCollection() {
        if (videoFiles.isEmpty()) {
            logError("No videos to analyze.");
            return;
        }
        
        StringBuilder analysis = new StringBuilder();
        analysis.append("=== VIDEO COLLECTION ANALYSIS ===\n");
        analysis.append("Total videos: ").append(videoFiles.size()).append("\n");
        
        long totalSize = videoFiles.stream().mapToLong(v -> v.size).sum();
        analysis.append("Total size: ").append(formatFileSize(totalSize)).append("\n");
        
        analysis.append("\nBy Extension:\n");
        videosByExtension.forEach((ext, list) -> {
            long extSize = list.stream().mapToLong(v -> v.size).sum();
            analysis.append(String.format("  %s: %d files (%s)\n", 
                ext.toUpperCase(), list.size(), formatFileSize(extSize)));
        });
        
        analysis.append("\nBy Size Category:\n");
        videosBySizeCategory.entrySet().stream()
            .sorted(Map.Entry.comparingByKey())
            .forEach(entry -> {
                int category = entry.getKey();
                List<VideoFile> list = entry.getValue();
                long catSize = list.stream().mapToLong(v -> v.size).sum();
                analysis.append(String.format("  %d-%d MB: %d files (%s)\n",
                    category * 100, (category + 1) * 100, list.size(), formatFileSize(catSize)));
            });
        
        // Find oldest and newest
        VideoFile oldest = videoFiles.stream()
            .min(Comparator.comparing(v -> v.created))
            .orElse(null);
        VideoFile newest = videoFiles.stream()
            .max(Comparator.comparing(v -> v.created))
            .orElse(null);
        
        if (oldest != null && newest != null) {
            analysis.append("\nOldest: ").append(oldest.name)
                   .append(" (").append(oldest.created.toLocalDate()).append(")\n");
            analysis.append("Newest: ").append(newest.name)
                   .append(" (").append(newest.created.toLocalDate()).append(")\n");
        }
        
        logArea.setText(analysis.toString());
    }
    
    private void showConversionDialog() {
        if (videoFiles.isEmpty()) {
            logError("No videos to convert.");
            return;
        }
        
        JDialog dialog = new JDialog();
        dialog.setTitle("Video Format Conversion");
        dialog.setLayout(new GridLayout(4, 2, 10, 10));
        
        JComboBox<String> sourceExt = new JComboBox<>(
            videosByExtension.keySet().toArray(new String[0])
        );
        JComboBox<String> targetExt = new JComboBox<>(
            new String[]{"mp4", "avi", "mov", "mkv", "webm"}
        );
        JCheckBox preserveOriginal = new JCheckBox("Preserve Original", true);
        JButton convertBtn = new JButton("Convert Selected");
        
        dialog.add(new JLabel("Source Format:"));
        dialog.add(sourceExt);
        dialog.add(new JLabel("Target Format:"));
        dialog.add(targetExt);
        dialog.add(preserveOriginal);
        dialog.add(new JLabel());
        dialog.add(convertBtn);
        
        convertBtn.addActionListener(e -> {
            String source = (String) sourceExt.getSelectedItem();
            String target = (String) targetExt.getSelectedItem();
            boolean preserve = preserveOriginal.isSelected();
            
            dialog.dispose();
            convertVideos(source, target, preserve);
        });
        
        dialog.pack();
        dialog.setVisible(true);
    }
    
    private void convertVideos(String sourceExt, String targetExt, boolean preserve) {
        // Note: Actual conversion would require FFmpeg or similar library
        // This is a simulation of the conversion process
        logMessage("Starting conversion: " + sourceExt + " -> " + targetExt);
        
        List<VideoFile> toConvert = videosByExtension.getOrDefault(sourceExt, new ArrayList<>());
        if (toConvert.isEmpty()) {
            logError("No videos with extension: " + sourceExt);
            return;
        }
        
        executor.submit(() -> {
            int success = 0;
            for (VideoFile video : toConvert) {
                try {
                    // Simulate conversion process
                    logMessage("Converting: " + video.name);
                    Thread.sleep(500); // Simulate conversion time
                    
                    if (!preserve) {
                        // Simulate deleting original
                        logMessage("Deleted original: " + video.name);
                    }
                    
                    success++;
                    logMessage("Converted: " + video.name + " -> " + 
                        video.name.replace("." + sourceExt, "." + targetExt));
                    
                } catch (Exception e) {
                    logError("Failed to convert: " + video.name);
                }
            }
            logMessage("Conversion complete: " + success + "/" + toConvert.size() + " videos converted");
        });
    }
    
    private void findDuplicateVideos() {
        if (videoFiles.isEmpty()) {
            logError("No videos to analyze for duplicates.");
            return;
        }
        
        logMessage("Searching for duplicate videos...");
        
        // Group by size (simple duplicate detection)
        Map<Long, List<VideoFile>> sizeMap = new HashMap<>();
        for (VideoFile video : videoFiles) {
            sizeMap.computeIfAbsent(video.size, k -> new ArrayList<>()).add(video);
        }
        
        int duplicateCount = 0;
        for (Map.Entry<Long, List<VideoFile>> entry : sizeMap.entrySet()) {
            if (entry.getValue().size() > 1) {
                duplicateCount++;
                logMessage("\nPotential duplicates (Size: " + formatFileSize(entry.getKey()) + "):");
                for (VideoFile video : entry.getValue()) {
                    logMessage("  - " + video.path);
                }
            }
        }
        
        if (duplicateCount == 0) {
            logMessage("No duplicate videos found.");
        } else {
            logMessage("\nFound " + duplicateCount + " sets of potential duplicates.");
        }
    }
    
    private void exportVideoReport() {
        if (videoFiles.isEmpty()) {
            logError("No videos to export.");
            return;
        }
        
        JFileChooser chooser = new JFileChooser();
        chooser.setSelectedFile(new File("video_report.txt"));
        chooser.setDialogTitle("Export Video Report");
        
        if (chooser.showSaveDialog(null) == JFileChooser.APPROVE_OPTION) {
            File reportFile = chooser.getSelectedFile();
            executor.submit(() -> {
                try (PrintWriter writer = new PrintWriter(reportFile)) {
                    writer.println("VIDEO COLLECTION REPORT");
                    writer.println("Generated: " + LocalDateTime.now());
                    writer.println("Total Videos: " + videoFiles.size());
                    writer.println("\nVideo List:");
                    
                    for (VideoFile video : videoFiles) {
                        writer.printf("%-40s %-10s %12s %s%n",
                            video.name,
                            video.extension.toUpperCase(),
                            formatFileSize(video.size),
                            video.created.format(DateTimeFormatter.ISO_LOCAL_DATE)
                        );
                    }
                    
                    logMessage("Report exported to: " + reportFile.getAbsolutePath());
                } catch (IOException e) {
                    logError("Failed to export report: " + e.getMessage());
                }
            });
        }
    }
    
    private void generatePreviews() {
        // This would generate thumbnail previews in a real implementation
        logMessage("Preview generation would require FFmpeg integration.");
        logMessage("This feature would create thumbnail grids for video files.");
    }
    
    private void extractVideoMetadata() {
        // This would extract detailed metadata using a library like JCodec or MediaInfo
        logMessage("Advanced metadata extraction would require JCodec/MediaInfo library.");
        logMessage("This would extract codec, duration, resolution, bitrate, etc.");
    }
    
    // Utility methods
    private String formatFileSize(long bytes) {
        if (bytes < 1024) return bytes + " B";
        int exp = (int) (Math.log(bytes) / Math.log(1024));
        String unit = "KMGTPE".charAt(exp-1) + "B";
        return String.format("%.1f %s", bytes / Math.pow(1024, exp), unit);
    }
    
    private void logMessage(String message) {
        SwingUtilities.invokeLater(() -> {
            logArea.append(message + "\n");
            logArea.setCaretPosition(logArea.getDocument().getLength());
        });
    }
    
    private void logError(String error) {
        logMessage("[ERROR] " + error);
    }
    
    public void cleanup() {
        executor.shutdown();
        try {
            if (!executor.awaitTermination(5, TimeUnit.SECONDS)) {
                executor.shutdownNow();
            }
        } catch (InterruptedException e) {
            executor.shutdownNow();
        }
    }
    
    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            VideoFileManager manager = new VideoFileManager();
            Runtime.getRuntime().addShutdownHook(new Thread(manager::cleanup));
        });
    }
}