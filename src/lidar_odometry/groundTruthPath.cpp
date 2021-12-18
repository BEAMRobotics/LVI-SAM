#include "utility.h"

class GroundTruthPath : public ParamServer
{
public:

    ros::Subscriber subGroundTruth;

    // map -> odom
    tf::Transform map_to_odom;
    tf::TransformBroadcaster tfMap2Odom;
    // odom -> base_link
    tf::TransformBroadcaster tfOdom2BaseLink;

    GroundTruthPath()
    {
        subGroundTruth = nh.subscribe<nav_msgs::Odometry>("/wamv/sensors/position/ground_truth_odometry", 5, &GroundTruthPath::groundTruthHandler, this, ros::TransportHints().tcpNoDelay());
        map_to_odom = tf::Transform(tf::createQuaternionFromRPY(0, 0, 0), tf::Vector3(0, 0, 0));
    }

    void groundTruthHandler(const nav_msgs::Odometry::ConstPtr& groundTruthMsg)
    {
        nav_msgs::Odometry odometry = groundTruthConverter(*groundTruthMsg);

        // publish transformation
        tf::Transform tCur;
        tf::poseMsgToTF(odometry.pose.pose, tCur);
        tf::StampedTransform odom_2_baselink = tf::StampedTransform(tCur, odometry.header.stamp, "odom", "base_link");
        tfOdom2BaseLink.sendTransform(odom_2_baselink);
    }
};


int main(int argc, char** argv)
{
    ros::init(argc, argv, "lidar");

    GroundTruthPath GTP;

    ROS_INFO("\033[1;32m----> Ground Truth Path Started.\033[0m");

    ros::spin();

    return 0;
}
