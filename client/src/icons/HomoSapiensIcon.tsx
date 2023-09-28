import HomoSapiens from './svg/homo_sapiens.svg';

interface IHomoSapiens {
  size: number;
}
function HomoSapiensIcon({ size }: IHomoSapiens) {
  return <img src={HomoSapiens} className="h-auto" style={{ width: size }} alt="homo sapiens icon" />;
}

export default HomoSapiensIcon;
