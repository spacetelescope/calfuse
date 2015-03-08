/*
 * GetTLE.java
 *
 * Created on January 16, 2005, 8:34 PM
 */
import java.net.URL;
import java.net.URLConnection;
import java.io.*;

import java.util.Map;
import java.util.List;
/**
 *
 * @author fred
 */
public class GetTLE {
    
    /** Creates a new instance of GetTLE */
    public GetTLE() {
    }
    
    /**
     * @param args the command line arguments
     */
    static final String host="www.space-track.org";
    static final String loginURL="/perl/login.pl";
    static final String dataURL="/perl/id_query.pl?ids=25791&timeframe=last5&common_name=yes&sort=catnum&descending=yes&ascii=yes&_submit=Submit&_submitted=1";
    static final String username="mromelfanger";
    static final String password="FuseJHU1";
  
    static final String outputFile = "/data1/fuse/calfuse/caltemp/five.tle";
    static final String logFile = "/data1/fuse/calfuse/caltemp/get_tle.logfile";
/*
    static final String outputFile = "/caltemp/five.tle";
    static final String logFile = "/caltemp/get_tle.logfile";
*/
    public static void main(String[] args) {
        // TODO code application logic here
        try {
            /*
             * Open the login url and get the session cookie back.
             */
            PrintWriter log = new PrintWriter(new OutputStreamWriter(new FileOutputStream(logFile)));
            PrintWriter out = new PrintWriter(new OutputStreamWriter(new FileOutputStream(outputFile)));
            
            URL u = new URL("http://"+host+loginURL+"?username="+username+"&password="+password+"&_submitted=1&_submit=Submit");

            URLConnection c = u.openConnection();
            Map m = c.getHeaderFields();
            List l = (List) m.get("Set-Cookie");
            String session = null;
            for(int i=0;i!=l.size();i++) {
                String s = (String) l.get(i);
                if(s.startsWith("spacetrack_session=")) {
                    int end = s.indexOf(";");
                    if(end < 0)
                        end = s.length();
                    session = s.substring(0, end);
                    break;
                }
            }
            InputStream is = c.getInputStream();
            BufferedReader br = new BufferedReader(new InputStreamReader(is));
            String s;
            while((s = br.readLine()) != null)
                log.println(s);
            is.close();
            
            if(session == null) {
                out.println("Could not get session id");
                System.exit(1);
            }

            /*
             * Use the session cookie to query for the fuse data.
             */
            
            u = new URL("http://"+host+dataURL);
            c = u.openConnection();
            c.setRequestProperty("Cookie", session);
            is = c.getInputStream();
            br = new BufferedReader(new InputStreamReader(is));
            /*
             * Read each line and keep what is between the pre's, add an extra
             * line at the end.
             */
            boolean collect = false;
	    out.println();
            while((s = br.readLine()) != null) {
                log.println(s);
                if(s.startsWith("<pre>")) {
/*                    br.readLine(); */
                    collect = true;
                } else if(s.startsWith("</pre>")) {
                    break;
                } else if(collect) {
                    out.println(s);
                }
            }
            out.println();
            br.close();
            is.close();
            out.close();
            log.close();

        } catch (Throwable e) {
            e.printStackTrace();
        }
    }
}
