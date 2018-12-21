/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package parser;

import java.io.File;

public class FileParseException extends Exception {
    public static FileParseException make(File file, String line, int lineIndex, String message) {
        return new FileParseException(file, line, lineIndex, message);
    }

    File file;
    String line;
    int lineIndex;

    public File getFile() {
        return file;
    }

    public String getLine() {
        return line;
    }

    public int getLineIndex() {
        return lineIndex;
    }

    /**
     * errorIndex == text.length は， あるべき行が存在しないことを表す．
     *
     * @param file      - file containing some error.
     * @param line      - line where the error exists
     * @param lineIndex - the row of the line in the file
     * @param message   - error message
     */
    FileParseException(File file, String line, int lineIndex, String message) {
        super(message);
        this.file = file;
        this.line = line;
        this.lineIndex = lineIndex;
    }

    @Override
    public String toString() {
        return "FileParseException{" +
                "file=" + file +
                ", line='" + line + '\'' +
                ", lineIndex=" + lineIndex +
                '}';
    }
}
