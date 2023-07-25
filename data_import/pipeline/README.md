




Lambda Configuration Steps
Pre-requisites: you have already set up an s3 bucket to store the h5ad files and created the Batch Job configuration

1. run `make build_lambda_layr` to create a lambda layer with the required python packages (.zip file)
2. create a new lambda layer in the aws console by uploading this .zip file
3. create a new lambda function in the aws console (Runtime: Python 3.10)
4. add the layer created in step 2 to the lambda function
5. copy the code from `lambda_function.py` into the lambda function
6. under "Configuration / Triggers" add a new trigger on the S3 Bucket for all object create events
7. Add a new Postgres Database Secret in the AWS Secrets Manager
8. Add the following policies to your lambda execution role

S3 Access, Log Creation, Batch Submit Job, Access to Secrets Manager(Inline Policy)
```json
{
    "Version": "2012-10-17",
    "Statement": [
        {
            "Effect": "Allow",
            "Action": [
                "logs:PutLogEvents",
                "logs:CreateLogGroup",
                "logs:CreateLogStream"
            ],
            "Resource": "arn:aws:logs:*:*:*"
        },
        {
            "Effect": "Allow",
            "Action": [
                "s3:GetObject"
            ],
            "Resource": "arn:aws:s3:::YOUR_S3_BUCKET/*"
        },
      {
            "Sid": "VisualEditor0",
            "Effect": "Allow",
            "Action": "batch:SubmitJob",
            "Resource": [
                "YOUR_JOB_DEFINITION_ARN:*",
                "YOUR_JOB_QUEUE_ARN"
            ]
        },
      {
    "Version": "2012-10-17",
    "Statement": [
        {
            "Sid": "VisualEditor0",
            "Effect": "Allow",
            "Action": [
                "secretsmanager:GetSecretValue",
                "secretsmanager:DescribeSecret"
            ],
            "Resource": "YOUR_SECRET_ARN"
        }
    ]
}
    ]
}
```

Also attach the following AWS provided policies: 
- AWSLambdaVPCAccessExecutionRole
- AWSBatchServiceEventTargetRole


9. Create a new security group for the lambda function with the following outbound rules:
- Type: Custom TCP, Protocol: TCP, Port Range: 5001, Destination: the security group of the EC2 instance
- Type: All Traffic, Protocol: All, Destination: 0.0.0.0/0

10. Add the security group created in step 9 to the EC2 instance as inbound rule source with protocol TCP and port range 5001
11. under "Configuration / VPC" add the VPC and subnets that the EC2 instance is in and add the security group created in step 9 to the lambda function


