import javax.imageio.ImageIO;
import javax.swing.*;
import java.awt.*;
import java.awt.image.BufferedImage;
import java.io.*;
import java.nio.file.*;
import java.util.*;
import java.util.List;
import java.util.concurrent.*;

public class ImageLoaderController {
    private List<File> imageFiles = new ArrayList<>();
    private int currentIndex = -1;
    private JLabel imageLabel;
    private JLabel infoLabel;
    private ExecutorService executor = Executors.newSingleThreadExecutor();
    
    // Supported image formats
    private static final Set<String> SUPPORTED_FORMATS = Set.of(
        "jpg", "jpeg", "png", "gif", "bmp", "webp"
    );
    
    public ImageLoaderController() {
        setupGUI();
    }
    
    private void setupGUI() {
        JFrame frame = new JFrame("Image Loader Controller");
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setLayout(new BorderLayout());
        
        // Control panel
        JPanel controlPanel = new JPanel();
        JButton loadBtn = new JButton("Load Images");
        JButton nextBtn = new JButton("Next");
        JButton prevBtn = new JButton("Previous");
        JButton zoomInBtn = new JButton("Zoom +");
        JButton zoomOutBtn = new JButton("Zoom -");
        
        controlPanel.add(loadBtn);
        controlPanel.add(prevBtn);
        controlPanel.add(nextBtn);
        controlPanel.add(zoomInBtn);
        controlPanel.add(zoomOutBtn);
        
        // Image display area
        imageLabel = new JLabel("No image loaded", SwingConstants.CENTER);
        imageLabel.setPreferredSize(new Dimension(600, 400));
        JScrollPane scrollPane = new JScrollPane(imageLabel);
        
        // Information panel
        infoLabel = new JLabel("Status: Ready");
        infoLabel.setBorder(BorderFactory.createEmptyBorder(5, 10, 5, 10));
        
        frame.add(controlPanel, BorderLayout.NORTH);
        frame.add(scrollPane, BorderLayout.CENTER);
        frame.add(infoLabel, BorderLayout.SOUTH);
        
        // Action listeners
        loadBtn.addActionListener(e -> loadImagesFromDirectory());
        nextBtn.addActionListener(e -> showNextImage());
        prevBtn.addActionListener(e -> showPreviousImage());
        zoomInBtn.addActionListener(e -> scaleImage(1.2));
        zoomOutBtn.addActionListener(e -> scaleImage(0.8));
        
        frame.pack();
        frame.setVisible(true);
    }
    
    private void loadImagesFromDirectory() {
        JFileChooser chooser = new JFileChooser();
        chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
        chooser.setDialogTitle("Select Image Directory");
        
        if (chooser.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) {
            File directory = chooser.getSelectedFile();
            executor.submit(() -> {
                try {
                    loadAllImages(directory);
                } catch (Exception e) {
                    SwingUtilities.invokeLater(() -> 
                        infoLabel.setText("Error: " + e.getMessage()));
                }
            });
        }
    }
    
    private void loadAllImages(File directory) throws IOException {
        imageFiles.clear();
        currentIndex = -1;
        
        if (!directory.exists() || !directory.isDirectory()) {
            throw new IOException("Invalid directory");
        }
        
        try (Stream<Path> paths = Files.walk(directory.toPath())) {
            paths.filter(Files::isRegularFile)
                .filter(this::isImageFile)
                .map(Path::toFile)
                .forEach(imageFiles::add);
        }
        
        SwingUtilities.invokeLater(() -> {
            if (imageFiles.isEmpty()) {
                infoLabel.setText("No images found in: " + directory.getName());
                imageLabel.setText("No images found");
            } else {
                infoLabel.setText("Loaded " + imageFiles.size() + " images");
                currentIndex = 0;
                loadAndDisplayImage(currentIndex);
            }
        });
    }
    
    private boolean isImageFile(Path path) {
        String fileName = path.getFileName().toString().toLowerCase();
        int dotIndex = fileName.lastIndexOf('.');
        if (dotIndex > 0) {
            String extension = fileName.substring(dotIndex + 1);
            return SUPPORTED_FORMATS.contains(extension);
        }
        return false;
    }
    
    private void loadAndDisplayImage(int index) {
        if (index < 0 || index >= imageFiles.size()) return;
        
        executor.submit(() -> {
            try {
                File imageFile = imageFiles.get(index);
                BufferedImage original = ImageIO.read(imageFile);
                
                if (original == null) {
                    SwingUtilities.invokeLater(() -> 
                        infoLabel.setText("Failed to load: " + imageFile.getName()));
                    return;
                }
                
                // Scale to fit display
                Image scaled = original.getScaledInstance(
                    Math.min(600, original.getWidth()),
                    Math.min(400, original.getHeight()),
                    Image.SCALE_SMOOTH
                );
                
                SwingUtilities.invokeLater(() -> {
                    imageLabel.setIcon(new ImageIcon(scaled));
                    imageLabel.setText(null);
                    infoLabel.setText(String.format(
                        "Image %d/%d: %s (%dx%d)",
                        index + 1, imageFiles.size(),
                        imageFile.getName(),
                        original.getWidth(), original.getHeight()
                    ));
                });
            } catch (IOException e) {
                SwingUtilities.invokeLater(() -> 
                    infoLabel.setText("Error loading: " + e.getMessage()));
            }
        });
    }
    
    private void showNextImage() {
        if (imageFiles.isEmpty()) return;
        currentIndex = (currentIndex + 1) % imageFiles.size();
        loadAndDisplayImage(currentIndex);
    }
    
    private void showPreviousImage() {
        if (imageFiles.isEmpty()) return;
        currentIndex = (currentIndex - 1 + imageFiles.size()) % imageFiles.size();
        loadAndDisplayImage(currentIndex);
    }
    
    private void scaleImage(double factor) {
        Icon icon = imageLabel.getIcon();
        if (icon instanceof ImageIcon) {
            ImageIcon imageIcon = (ImageIcon) icon;
            Image original = imageIcon.getImage();
            int newWidth = (int)(original.getWidth(null) * factor);
            int newHeight = (int)(original.getHeight(null) * factor);
            
            Image scaled = original.getScaledInstance(
                Math.max(50, newWidth),
                Math.max(50, newHeight),
                Image.SCALE_SMOOTH
            );
            
            imageLabel.setIcon(new ImageIcon(scaled));
            infoLabel.setText(String.format("Scaled to: %dx%d", newWidth, newHeight));
        }
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
            ImageLoaderController controller = new ImageLoaderController();
            Runtime.getRuntime().addShutdownHook(new Thread(controller::cleanup));
        });
    }
}