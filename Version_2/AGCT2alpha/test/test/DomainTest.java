package test;

import junit.framework.Assert;
import org.junit.Before;
import org.junit.Test;
import parser.FileParseException;

import java.io.File;
import java.io.FileNotFoundException;

public class DomainTest {

    @Before
    public void init() {
        AGCT.initAGCT();
    }


    @Test
    public void testReadFile() {
        File file = new File("AGCT2alpha/sample/raw/__Test_570.txt");
        try {
            Domain.getInstance().readFile(file);
        } catch (FileParseException e) {
            Assert.fail();
        } catch (FileNotFoundException e) {
            Assert.fail();
        }
    }

    @Test
    public void testReadFileWrongSpecies() {
        readFileExpectingParseError("Species");
    }

    @Test
    public void testReadFileWrongEntity() {
        readFileExpectingParseError("Entity");
    }

    @Test
    public void testReadFileWrongEntityGroup() {
        readFileExpectingParseError("EntityGroup");
    }

    @Test
    public void testReadFileWrongReplicate() {
        readFileExpectingParseError("Species");
    }

    @Test
    public void testReadFileWrongTime() {
        readFileExpectingParseError("Species");
    }

    @Test
    public void testReadFileWrongData() {
        readFileExpectingParseError("Data");
    }

    private void readFileExpectingParseError(String cause) {
        File file = new File("AGCT2alpha/sample/raw/__Test_570_wrong" + cause + ".txt");
        try {
            Domain.getInstance().readFile(file);
        } catch (FileParseException e) {
            return;
        } catch (FileNotFoundException e) {
            Assert.fail();
        }
        Assert.fail();
    }

    @Test
    public void testReadFileNoFile() {
        File file = new File("AGCT2alpha/sample/raw/no_such_file");
        try {
            Domain.getInstance().readFile(file);
        } catch (FileParseException e) {
            Assert.fail();
        } catch (FileNotFoundException e) {
            return;
        }
        Assert.fail();
    }
}
