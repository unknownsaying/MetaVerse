import java.io.*;
import java.nio.charset.StandardCharsets;

public class IndeedIOStream implements AutoCloseable {
    private final DataInputStream input;
    private final DataOutputStream output;
    
    public IndeedIOStream(InputStream in, OutputStream out) {
        this.input = new DataInputStream(in);
        this.output = new DataOutputStream(out);
    }
    
    // Reading methods
    public int readInt() throws IOException {
        return input.readInt();
    }
    
    public long readLong() throws IOException {
        return input.readLong();
    }
    
    public double readDouble() throws IOException {
        return input.readDouble();
    }
    
    public boolean readBoolean() throws IOException {
        return input.readBoolean();
    }
    
    public String readUTF() throws IOException {
        return input.readUTF();
    }
    
    public byte[] readBytes(int length) throws IOException {
        byte[] buffer = new byte[length];
        input.readFully(buffer);
        return buffer;
    }
    
    // Writing methods
    public void writeInt(int value) throws IOException {
        output.writeInt(value);
        output.flush();
    }
    
    public void writeLong(long value) throws IOException {
        output.writeLong(value);
        output.flush();
    }
    
    public void writeDouble(double value) throws IOException {
        output.writeDouble(value);
        output.flush();
    }
    
    public void writeBoolean(boolean value) throws IOException {
        output.writeBoolean(value);
        output.flush();
    }
    
    public void writeUTF(String value) throws IOException {
        output.writeUTF(value);
        output.flush();
    }
    
    public void writeBytes(byte[] data) throws IOException {
        output.write(data);
        output.flush();
    }
    
    // Convenience methods
    public void writeAscii(String text) throws IOException {
        byte[] bytes = text.getBytes(StandardCharsets.US_ASCII);
        writeInt(bytes.length);
        writeBytes(bytes);
    }
    
    public String readAscii() throws IOException {
        int length = readInt();
        byte[] bytes = readBytes(length);
        return new String(bytes, StandardCharsets.US_ASCII);
    }
    
    // Buffer management
    public void flush() throws IOException {
        output.flush();
    }
    
    public int available() throws IOException {
        return input.available();
    }
    
    // Resource cleanup
    @Override
    public void close() throws IOException {
        try {
            input.close();
        } finally {
            output.close();
        }
    }
    
    // Demo usage
    public static void main(String[] args) {
        try {
            ByteArrayOutputStream byteOut = new ByteArrayOutputStream();
            IndeedIOStream stream = new IndeedIOStream(
                new ByteArrayInputStream(byteOut.toByteArray()), 
                byteOut
            );
            
            // Write test data
            stream.writeInt(42);
            stream.writeUTF("Hello Indeed!");
            stream.writeBoolean(true);
            
            // Read back
            stream = new IndeedIOStream(
                new ByteArrayInputStream(byteOut.toByteArray()),
                null
            );
            
            System.out.println("Int: " + stream.readInt());
            System.out.println("String: " + stream.readUTF());
            System.out.println("Boolean: " + stream.readBoolean());
            
            stream.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}
// Real-world usage example
try (IndeedIOStream io = new IndeedIOStream(socket.getInputStream(), socket.getOutputStream())) {
    io.writeUTF("Client Request");
    String response = io.readUTF();
    // Process response...
}